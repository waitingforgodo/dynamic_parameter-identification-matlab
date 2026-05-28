function generate_load_regressor(path_to_urdf)
% ----------------------------------------------------------------------
% The function generates regressor for the load, which is assumed to be
% a rigid body
% -----------------------------------------------------------------------
% Parse urdf to get robot description
ur10 = parse_urdf(path_to_urdf);

% Create symbolic generilized coordiates, their first and second deriatives
q_sym = sym('q%d',[7,1],'real');
qd_sym = sym('qd%d',[7,1],'real');
q2d_sym = sym('q2d%d',[7,1],'real');

% ------------------------------------------------------------------------
% Getting gradient of energy functions, to derive dynamics
% ------------------------------------------------------------------------
T_pk = sym(zeros(4,4,7)); % transformation between links
w_kk(:,1) = sym(zeros(3,1)); % angular velocity k in frame k
v_kk(:,1) = sym(zeros(3,1)); % linear velocity of the origin of frame k in frame k
g_kk(:,1) = sym([0,0,9.81])'; % vector of graviatational accelerations in frame k
p_kk(:,1) = sym(zeros(3,1)); % origin of frame k in frame k

for i = 1:7
    jnt_axs_k = str2num(ur10.robot.joint{i}.axis.Attributes.xyz)';
    % Transformation from parent link frame p to current joint frame
    rpy_k = sym(str2num(ur10.robot.joint{i}.origin.Attributes.rpy));
    R_pj = RPY(rpy_k);
    R_pj(abs(R_pj)<sqrt(eps)) = sym(0); % to avoid numerical errors
    p_pj = str2num(ur10.robot.joint{i}.origin.Attributes.xyz)';
    T_pj = sym([R_pj, p_pj; zeros(1,3), 1]); % to avoid numerical errors
    % Tranformation from joint frame of the joint that rotaties body k to
    % link frame. The transformation is pure rotation
    R_jk = Rot(q_sym(i),sym(jnt_axs_k));
    p_jk = sym(zeros(3,1));
    T_jk = [R_jk, p_jk; sym(zeros(1,3)),sym(1)];
    % Transformation from parent link frame p to current link frame k
    T_pk(:,:,i) = T_pj*T_jk;
    z_kk(:,i) = sym(jnt_axs_k);
        
    w_kk(:,i+1) = T_pk(1:3,1:3,i)'*w_kk(:,i) + sym(jnt_axs_k)*qd_sym(i);
    v_kk(:,i+1) = T_pk(1:3,1:3,i)'*(v_kk(:,i) + cross(w_kk(:,i),sym(p_pj)));
    g_kk(:,i+1) = T_pk(1:3,1:3,i)'*g_kk(:,i);
    p_kk(:,i+1) = T_pk(1:3,1:3,i)'*(p_kk(:,i) + sym(p_pj));
        
    beta_K(i,:) = [sym(0.5)*w2wtlda(w_kk(:,i+1)),...
                   v_kk(:,i+1)'*vec2skewSymMat(w_kk(:,i+1)),...
                   sym(0.5)*v_kk(:,i+1)'*v_kk(:,i+1)];
    beta_P(i,:) = [sym(zeros(1,6)), g_kk(:,i+1)',...
                   g_kk(:,i+1)'*p_kk(:,i+1)];
end

% --------------------------------------------------------------------
% Gradient of the kinetic and potential energy of the load
% --------------------------------------------------------------------
% Transformation from link 7 frame to end-effector frame
rpy_ee = sym(str2num(ur10.robot.joint{8}.origin.Attributes.rpy));
R_7ee = RPY(rpy_ee);
R_7ee(abs(R_7ee)<sqrt(eps)) = sym(0); % to avoid numerical errors
p_7ee = str2num(ur10.robot.joint{8}.origin.Attributes.xyz)';
T_7ee = sym([R_7ee, p_7ee; zeros(1,3), 1]); % to avoid numerical errors

w_eeee = T_7ee(1:3,1:3)'*w_kk(:,8);
v_eeee = T_7ee(1:3,1:3)'*(v_kk(:,8) + cross(w_kk(:,i+1),sym(p_7ee)));
g_eeee = T_7ee(1:3,1:3)'*g_kk(:,8);
p_eeee = T_7ee(1:3,1:3)'*(p_kk(:,8) + sym(p_7ee));

beta_Kl = [sym(0.5)*w2wtlda(w_eeee), v_eeee'*vec2skewSymMat(w_eeee),...
            sym(0.5)*(v_eeee'*v_eeee)];
        
beta_Pl = [sym(zeros(1,7)), g_eeee', g_eeee'*p_eeee];


% ---------------------------------------------------------------------
% Dynamic regressor of the load
% ---------------------------------------------------------------------
beta_Ll = beta_Kl - beta_Pl;
dbetaLl_dq = jacobian(beta_Ll,q_sym)';
dbetaLl_dqd = jacobian(beta_Ll,qd_sym)';
tl = sym(zeros(7,10));
for i = 1:7
   tl = tl + diff(dbetaLl_dqd,q_sym(i))*qd_sym(i)+...
                diff(dbetaLl_dqd,qd_sym(i))*q2d_sym(i);
end
Y_l = tl - dbetaLl_dq;

% Generate a function from a symbolic expression
matlabFunction(Y_l,'File','autogen/load_regressor_UR10E',...
               'Vars',{q_sym, qd_sym, q2d_sym});
