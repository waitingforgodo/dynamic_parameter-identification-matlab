function TestBaseParams()
clc
% ---------------------------------------------------------------------
% Test if base parameters are found correctly by comparing
% torque prediction of the model with standard parameters with 
% the torque prediction of the model with base parameters
% ----------------------------------------------------------------------
addpath(genpath('D:\myWork')); 
robot = importrobot('Marvin.urdf');

joint = 7;
Marvin_pi = zeros(70, joint);
for i = 1:joint
currentBody = robot.Bodies{i};
    m = currentBody.Mass;
    com = currentBody.CenterOfMass;
    inertiaVec = currentBody.Inertia;
    localParams = [m; com(:); inertiaVec(:)];
    startIdx = (i - 1) * 10 + 1;
    endIdx = i * 10;
    Marvin_pi(startIdx:endIdx) = localParams;
end

model_type = 'Newton-Euler';
% model_type = 'Lagrangian';

[~, baseQR] = GenBaseParameters(robot,model_type);

bb = baseQR.numberOfBaseParameters;
E = baseQR.permutationMatrix;
beta = baseQR.beta;

% Position, velocity and acceleration limits
q_min = -pi*ones(joint,1);
q_max = pi*ones(joint,1);
qd_max = 3*pi*ones(joint,1);
q2d_max = 6*pi*ones(joint,1);

% On random positions, velocities, aceeleations
for i = 1:100
    q_rnd = q_min + (q_max - q_min).*rand(joint,1);
    qd_rnd = -qd_max + 2*qd_max.*rand(joint,1);
    q2d_rnd = -q2d_max + 2*q2d_max.*rand(joint,1);
    
    if strcmpi(model_type, 'Newton-Euler')
        Yi = SymForm_NT_Y(q_rnd,qd_rnd,q2d_rnd);
    elseif strcmpi(model_type, 'Lagrangian')
        Yi = SymForm_LG_Y(q_rnd,qd_rnd,q2d_rnd);
    else
        assert(false,'Model type error');
    end

    tau_full = Yi*Marvin_pi;
    
    pi_lgr_base = [eye(bb) beta]*E'*Marvin_pi;
    Y_base = Yi*E(:,1:bb);
    tau_base = Y_base*pi_lgr_base;
    assert(norm(tau_full - tau_base) < 1e-6);
end

msg = model_type + " Model Rigid Body Base Dynamics Test - OK!";
fprintf(msg);