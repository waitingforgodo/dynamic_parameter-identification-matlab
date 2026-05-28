% clear all; close all; clc;

% Define path to a urdf file
path_to_urdf = 'Marvin.urdf';

% 确保并行池已启动
if isempty(gcp('nocreate'))
    parpool;
end
% Generate functions for dynamics based on Lagrange method
% Note that it might take some time
% generate_rb_dynamics();
% generate_friction_eq();


% Generate regressors for inverse dynamics of the robot, friction and load
% Note that it might take some time
% generate_rb_regressor();
%generate_load_regressor(path_to_urdf);


% Run tests
% test_rb_inverse_dynamics()
% test_base_params()


% Perform QR decompostion in order to get base parameters of the robot
include_motor_dynamics = 0;
cacheFile = 'baseQR_7dof.mat';
[pi_lgr_base, baseQR] = base_params_qr(include_motor_dynamics, cacheFile);


% % Estimate drive gains
% % drive_gains = estimate_drive_gains(baseQR, 'PC-OLS');
% % Or use those found in the paper by De Luca
% % drive_gains = [14.87; 13.26; 11.13; 10.62; 11.03; 11.47]; 
% 
% 
% % Estimate dynamic parameters
% path_to_est_data = '2026_05_26_120_5_torque.csv';      idxs = [2, 110];
% 
% % 数据截取窗口 [t_start, t_end]（秒）。T=120 轨迹建议跳过前 2s 启动段。
% % path_to_est_data = out;  idxs = [2, 110];
% 
% % idx_mode = 'time';
% % sol = estimate_dynamic_params(path_to_est_data, idxs, ...
% %                               baseQR, 'URDF-REFINE', [], idx_mode);
% % sol = estimate_dynamic_params(path_to_est_data, idxs, ...
% %                               baseQR, 'PC-OLS', [], idx_mode);
% sol = estimate_dynamic_params(path_to_est_data, idxs, ...
%                               baseQR, 'OLS', [], idx_mode);
% 
% % Validate estimated parameters（验证段可与辨识段相同或错开）
% % path_to_val_data = out; idxs = [2, 110];
% 
% path_to_val_data='2026_05_26_120_5_torque.csv';  idxs = [2, 110];
% 
% rre = validate_dynamic_params(path_to_val_data, idxs, ...
%                                baseQR, sol.pi_b, sol.pi_fr, idx_mode);
% 
% disp("辨识结果rre值:")
% disp(rre)










