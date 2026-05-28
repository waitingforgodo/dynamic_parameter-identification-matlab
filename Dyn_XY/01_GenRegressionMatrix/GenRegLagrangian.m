function GenRegLagrangian(robot)
%% Lagrangian公式推导回归矩阵 
% tau = Y(q, dq, ddq) pi

    %% 符号  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nJoints = 7;
    q_sym = sym('q%d',[nJoints,1],'real');
    dq_sym = sym('dq%d',[nJoints,1],'real');
    ddq_sym = sym('ddq%d',[nJoints,1],'real');

    T_pk = sym(zeros(4,4,nJoints)); % transformation between links
    w_kk(:,1) = sym(zeros(3,1)); % angular velocity k in frame k
    v_kk(:,1) = sym(zeros(3,1)); % linear velocity of the origin of frame k in frame k
    g_kk(:,1) = sym([0,0,9.81])'; % vector of graviatational accelerations in frame k
    p_kk(:,1) = sym(zeros(3,1)); % origin of frame k in frame k
    Z = [0 0 1.0]';

     %% 计算  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 运动学Fk  
    [T_pk(:,:,1),T_pk(:,:,2),T_pk(:,:,3),T_pk(:,:,4),T_pk(:,:,5),T_pk(:,:,6),T_pk(:,:,7)] = RobotFk(robot,q_sym(1),q_sym(2),q_sym(3),q_sym(4),q_sym(5),q_sym(6),q_sym(7));  

    for i = 1:nJoints   
        p_pj = T_pk(1:3,4,i);
        if i < nJoints
            w_kk(:,i+1) = T_pk(1:3,1:3,i)'*w_kk(:,i) + sym(Z)*dq_sym(i);
            v_kk(:,i+1) = T_pk(1:3,1:3,i)'*(v_kk(:,i) + cross(w_kk(:,i),sym(p_pj)));
        else
            w_kk(:,i+1) = T_pk(1:3,1:3,i)'*w_kk(:,i);
            v_kk(:,i+1) = T_pk(1:3,1:3,i)'*(v_kk(:,i) + cross(w_kk(:,i),sym(p_pj))) + sym(Z)*dq_sym(i);
        end
    
        g_kk(:,i+1) = T_pk(1:3,1:3,i)'*g_kk(:,i);
        p_kk(:,i+1) = T_pk(1:3,1:3,i)'*(p_kk(:,i) + sym(p_pj));
    
        beta_K(i,:) = [sym(0.5)*w2wtlda(w_kk(:,i+1)),...
                       v_kk(:,i+1)'*vec2skewSymMat(w_kk(:,i+1)),...
                       sym(0.5)*v_kk(:,i+1)'*v_kk(:,i+1)];
        beta_P(i,:) = [sym(zeros(1,6)), g_kk(:,i+1)',...
                       g_kk(:,i+1)'*p_kk(:,i+1)];
    end
    
    beta_Lf = [beta_K(1,:)-beta_P(1,:), beta_K(2,:)-beta_P(2,:), beta_K(3,:)-beta_P(3,:), beta_K(4,:)-beta_P(4,:),...
               beta_K(5,:)-beta_P(5,:), beta_K(6,:)-beta_P(6,:), beta_K(7,:)-beta_P(7,:)];

    dbetaLf_dq = jacobian(beta_Lf,q_sym)';
    dbetaLf_dqd = jacobian(beta_Lf,dq_sym)';

    tf = sym(zeros(nJoints,nJoints*10));
    for i = 1:nJoints
       tf = tf + diff(dbetaLf_dqd,q_sym(i))*dq_sym(i) + diff(dbetaLf_dqd,dq_sym(i))*ddq_sym(i);
    end

    Y = tf - dbetaLf_dq;
    
    % F = Friction(dq_sym);

    % Generate function from a symbolic expression for the regressor
    matlabFunction(Y,'File','D:\myWork\Dyn_XY\GenRegressionMatrix\SymbolicFormula\SymForm_LG_Y',...
                   'Vars',{q_sym,dq_sym,ddq_sym});

    % matlabFunction(F,'File','D:\SymForm_F','Vars',{dq_sym});

end

%% 三维向量转反对称矩阵
function m = vec2skewSymMat(v)
    m = [    0,  -v(3),  v(2);
          v(3),      0, -v(1);
         -v(2),   v(1),    0];
end

%% 三维向量映射成包含其二阶项的六维行向量
function out = w2wtlda(w)
    out = [w(1)^2, 2*w(1)*w(2), 2*w(1)*w(3),...
            w(2)^2, 2*w(2)*w(3), w(3)^2];
end

%% 摩擦力 
function [frctn] = Friction(dq)
% ----------------------------------------------------------------------
% The function computes friction regressor for each joint of the robot.
% Fv*qd + Fc*sign(qd) + F0
% ---------------------------------------------------------------------
    nJoints = 7;
    %frctn = zeros(nJoints, nJoints*3);
    for i = 1:nJoints
        frctn(i,3*i-2:3*i) = [dq(i), sign(dq(i)), 1];
    end
end