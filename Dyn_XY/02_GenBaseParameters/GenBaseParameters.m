function [pi_lgr_base, baseQR] = GenBaseParameters(robot,model_type)

clear sym

% ----------------------------------------------------------------------
% Seed the random number generator based on the current time
rng('shuffle');


% ------------------------------------------------------------------------
% Set limits on posistion and velocities
% ------------------------------------------------------------------------
num = 7;
q_min = [robot.Bodies{1,1}.Joint.PositionLimits(1) robot.Bodies{1,2}.Joint.PositionLimits(1) robot.Bodies{1,3}.Joint.PositionLimits(1) ...
         robot.Bodies{1,4}.Joint.PositionLimits(1) robot.Bodies{1,5}.Joint.PositionLimits(1) robot.Bodies{1,6}.Joint.PositionLimits(1) ...
         robot.Bodies{1,7}.Joint.PositionLimits(1)]';
q_max = [robot.Bodies{1,1}.Joint.PositionLimits(2) robot.Bodies{1,2}.Joint.PositionLimits(2) robot.Bodies{1,3}.Joint.PositionLimits(2) ...
         robot.Bodies{1,4}.Joint.PositionLimits(2) robot.Bodies{1,5}.Joint.PositionLimits(2) robot.Bodies{1,6}.Joint.PositionLimits(2) ...
         robot.Bodies{1,7}.Joint.PositionLimits(2)]';
dq_max = 3*q_max;
ddq_max = 6*q_max;

% -----------------------------------------------------------------------
% Standard dynamics paramters of the robot in symbolic form
% -----------------------------------------------------------------------
im  = sym('im%d',[num,1],'real');
m   = sym('m%d', [num,1],'real');
hx  = sym('px%d',[num,1],'real');
hy  = sym('py%d',[num,1],'real');
hz  = sym('pz%d',[num,1],'real');
ixx = sym('Ixx%d',[num,1],'real');
ixy = sym('Ixy%d',[num,1],'real');
ixz = sym('Ixz%d',[num,1],'real');
iyy = sym('Iyy%d',[num,1],'real');
iyz = sym('Iyz%d',[num,1],'real');
izz = sym('Izz%d',[num,1],'real');

% Vector of symbolic parameters
for i = 1:num
    pi_lgr_sym(:,i) = [m(i),hx(i),hy(i),hz(i),ixx(i),ixy(i),ixz(i),iyy(i),iyz(i),izz(i)]';
end
[nLnkPrms, nLnks] = size(pi_lgr_sym);
pi_lgr_sym = reshape(pi_lgr_sym, [nLnkPrms*nLnks, 1]);


% -----------------------------------------------------------------------
% Find relation between independent columns and dependent columns
% -----------------------------------------------------------------------
% Get observation matrix of identifiable paramters
nSample = 50;
W = [];
for i = 1:nSample
    q_rnd = q_min + (q_max - q_min).*rand(num,1);
    dq_rnd = -dq_max + 2*dq_max.*rand(num,1);
    ddq_rnd = -ddq_max + 2*ddq_max.*rand(num,1);
    
    if strcmpi(model_type, 'Newton-Euler')
        Y = SymForm_NT_Y(q_rnd, dq_rnd, ddq_rnd);
    elseif strcmpi(model_type, 'Lagrangian')
        Y = SymForm_LG_Y(q_rnd, dq_rnd, ddq_rnd);
    else
        assert(false,'Model type error');
    end

    W = vertcat(W,Y);
end

% QR decomposition with pivoting: W*E = Q*R
%   R is upper triangular matrix
%   Q is unitary matrix
%   E is permutation matrix
[Q, R, E] = qr(W);

% matrix W has rank bb which is number of base parameters 
bb = rank(W);

% R = [R1 R2; 
%      0  0]
% R1 is bbxbb upper triangular and reguar matrix
% R2 is bbx(c-bb) matrix where c is number of standard parameters
R1 = R(1:bb,1:bb);
R2 = R(1:bb,bb+1:end);
beta = R1\R2; % the zero rows of K correspond to independent columns of WP
beta(abs(beta)<sqrt(eps)) = 0; % get rid of numerical errors
% W2 = W1*beta

% Make sure that the relation holds
W1 = W*E(:,1:bb);
W2 = W*E(:,bb+1:end);

residual = norm(W2 - W1*beta);
relative_error = residual / (norm(W2) + eps);
% LG 表达式更长，W 浮点累积误差通常比 NE 大，1e-9 过严
relTol = 1e-8;
assert(relative_error < relTol, ...
    sprintf('Found relationship between W1 and W2 is not correct. Relative error: %e (tol=%e)', ...
    relative_error, relTol));
% assert(norm(W2 - W1*beta) < 1e-6,... 
        % 'Found realationship between W1 and W2 is not correct\n');

% -----------------------------------------------------------------------
% Find base parmaters
% -----------------------------------------------------------------------
pi1 = E(:,1:bb)'*pi_lgr_sym; % independent paramters
pi2 = E(:,bb+1:end)'*pi_lgr_sym; % dependent paramteres

% all of the expressions below are equivalent
pi_lgr_base = pi1 + beta*pi2;
% pi_lgr_base = [eye(bb) beta]*[pi1;pi2];
% pi_lgr_base = [eye(bb) beta]*E'*pi_lgr_sym;

% ---------------------------------------------------------------------
% Create structure with the result of QR decompositon a
% ---------------------------------------------------------------------
baseQR = struct;
baseQR.numberOfBaseParameters = bb;
baseQR.permutationMatrix = E;
baseQR.beta = beta;
baseQR.motorDynamicsIncluded = 0;
