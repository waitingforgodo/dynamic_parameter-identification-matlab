function generate_friction_eq()
% 每个关节需要 3 个参数：
% 粘性摩擦系数 B
% 库伦摩擦系数 Fc
% 静摩擦 / 偏置摩擦 Fs
% tau_frcn=B.dot{q}_{粘性摩擦}
% + F_c.{sign}(dot{q})_{{库伦摩擦}}
% + F_s{静摩擦} tau_frcn 是 7×1 向量，对应7个关节的摩擦力矩。
% Create symbolic generilized coordiates, their first and second deriatives
qd_sym = sym('qd%d',[7,1],'real');

% Create symbolic experssions for friction parameters
pi_frcn = sym('pi_frcn_%d%d', [21,1], 'real');

% Friction torque
% pi_frcn_tmp = reshape(pi_frcn, [3, 7])';
pi_frcn_tmp = reshape(pi_frcn, 7, 3);  % 直接 7行3列，绝对安全
tau_frcn = pi_frcn_tmp(:,1).*qd_sym + pi_frcn_tmp(:,2).*sign(qd_sym) + pi_frcn_tmp(:,3);

% Generate a fucnction from symbolic expressions
matlabFunction(tau_frcn, 'File','autogen/F_vctr_fcn',...
               'Vars',{qd_sym, pi_frcn}, 'Optimize', true);