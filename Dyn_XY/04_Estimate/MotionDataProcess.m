function [md] = MotionDataProcess()

%% 拉取原始数据
motionData = readtable('D:\myWork\Dyn_XY\Scope.csv', 'VariableNamingRule', 'preserve');
[nRow, ~] = size(motionData);
Pos = motionData{2:nRow, 1:7};
Vel = motionData{2:nRow, 8:14};
Tor = motionData{2:nRow, 22:28};
[nRow, ~] = size(Pos);

%% 通用参数
dt = 0.001;           % 采样周期 (s)
trim_edge_s = 1.0;    % 两端裁剪时长 (s)
nTrimEdge = round(trim_edge_s / dt);

%% 速度滤波（辨识直接使用 VelFlt）
VelFlt = ProcessVelFilter(Vel, dt);

%% 转矩滤波（辨识直接使用 TorFlt）
TorFlt = ProcessTorFilter(Tor, dt);

%% 加速度处理
Acc_final = ProcessAccFromPos(Pos, dt);

%% 对齐与两端裁剪
[idxUse, tCrop, Pos_crop, Vel_crop, VelFlt_crop, Acc_crop, Tor_crop, TorFlt_crop] = ...
    alignAndCropMotionData(Pos, Vel, VelFlt, Acc_final, Tor, TorFlt, dt, nTrimEdge);

fprintf('辨识数据: n=%d, t=%.3f~%.3fs\n', numel(idxUse), tCrop(1), tCrop(end));

%% 打印输出
plotMotionData(tCrop, Pos_crop, Vel_crop, VelFlt_crop, Acc_crop, Tor_crop, TorFlt_crop);

%% 输出辨识结构体
deg2rad = pi / 180;

md = struct();
md.q = Pos_crop * deg2rad;
md.dq = VelFlt_crop * deg2rad;
md.ddq = Acc_crop * deg2rad;
md.tor = TorFlt_crop;

end

%% ========== 加速度处理 ==========

function Acc = ProcessAccFromPos(Pos, dt, opts)
% 从位置数据计算各关节加速度
% 流程：Hampel → 强低通+SG → 中心差分 → 速度轻平滑 → 中心差分 → 尖峰修复

if nargin < 3, opts = struct(); end

cfg.order_strong = 3;
cfg.order_vel = 3;
cfg.hampel_k = 5;
cfg.hampel_nsigma = 3;
cfg.acc_trim_start = 100;
cfg.fc_strong_j         = [5, 5, 5, 5, 5, 3, 3];
cfg.framelen_strong_j   = [51, 51, 51, 51, 51, 101, 101];
cfg.framelen_vel_j      = [11, 11, 11, 11, 11, 31, 31];
cfg.acc_ref_fc_j        = [3, 3, 3, 3, 3, 2, 2];
cfg.acc_repair_nsigma_j = [2.0, 2.0, 2.0, 2.0, 2.0, 1.3, 1.3];
fn = fieldnames(opts);
for k = 1:numel(fn), cfg.(fn{k}) = opts.(fn{k}); end
cfg.framelen_strong_j = ensureOdd(cfg.framelen_strong_j);
cfg.framelen_vel_j = ensureOdd(cfg.framelen_vel_j);

[nRow, nJoint] = size(Pos);
fs = 1 / dt;
nAcc = nRow - 4;
Acc = zeros(nAcc, nJoint);

for i = 1:nJoint
    fc_strong = cfg.fc_strong_j(i);
    framelen_strong = cfg.framelen_strong_j(i);
    framelen_vel = cfg.framelen_vel_j(i);

    pos_filt = designfilt('lowpassiir', 'FilterOrder', 4, ...
        'HalfPowerFrequency', fc_strong / (fs / 2), 'DesignMethod', 'butter');

    pos_clean = hampel(Pos(:, i), cfg.hampel_k, cfg.hampel_nsigma);
    pos_lp = filtfilt(pos_filt, pos_clean);
    pos_strong = sgolayfilt(pos_lp, cfg.order_strong, framelen_strong);

    vel_pos = (pos_strong(3:end) - pos_strong(1:end-2)) / (2 * dt);
    vel_light = sgolayfilt(vel_pos, cfg.order_vel, framelen_vel);

    acc_raw = (vel_light(3:end) - vel_light(1:end-2)) / (2 * dt);
    acc_h = hampel(acc_raw, cfg.hampel_k, cfg.hampel_nsigma);
    Acc(:, i) = repairAccSpikes(acc_h, dt, cfg.acc_ref_fc_j(i), ...
        cfg.acc_repair_nsigma_j(i), cfg.acc_trim_start);
end

end

function v = ensureOdd(v)
    v(mod(v, 2) == 0) = v(mod(v, 2) == 0) - 1;
end

function accOut = repairAccSpikes(accIn, dt, fcRef, nsigma, nTrimStart)
    fs = 1 / dt;
    lp = designfilt('lowpassiir', 'FilterOrder', 2, ...
        'HalfPowerFrequency', fcRef / (fs / 2), 'DesignMethod', 'butter');
    accRef = filtfilt(lp, accIn);
    dev = accIn - accRef;
    scale = max(1.4826 * mad(dev, 1), eps);
    mask = abs(dev) > nsigma * scale;
    accOut = accIn;
    accOut(mask) = accRef(mask);
    if nTrimStart > 0
        n = min(nTrimStart, numel(accOut));
        accOut(1:n) = accRef(1:n);
    end
end

%% ========== 速度滤波 ==========

function VelFlt = ProcessVelFilter(Vel, dt, opts)
% 速度滤波：Hampel 清野值 → Butterworth 低通（带宽与 ddq 对齐）

if nargin < 3, opts = struct(); end

cfg.filter_order = 6;
cfg.hampel_k = 5;
cfg.hampel_nsigma = 3;
cfg.fc_vel_j = [4, 4, 4, 4, 4, 3, 3];
fn = fieldnames(opts);
for k = 1:numel(fn), cfg.(fn{k}) = opts.(fn{k}); end

[nRow, nJoint] = size(Vel);
fs = 1 / dt;
VelFlt = zeros(nRow, nJoint);

for i = 1:nJoint
    fc = cfg.fc_vel_j(i);
    vel_filt = designfilt('lowpassiir', 'FilterOrder', cfg.filter_order, ...
        'HalfPowerFrequency', fc / (fs / 2), 'DesignMethod', 'butter');

    vel_clean = hampel(Vel(:, i), cfg.hampel_k, cfg.hampel_nsigma);
    VelFlt(:, i) = filtfilt(vel_filt, vel_clean);
end

end

%% ========== 转矩滤波 ==========

function TorFlt = ProcessTorFilter(Tor, dt, opts)
% 转矩滤波：Hampel 清野值 → Butterworth 低通 → 可选 SG 轻平滑

if nargin < 3, opts = struct(); end

cfg.filter_order = 6;
cfg.order_sg = 3;
cfg.hampel_k = 5;
cfg.hampel_nsigma = 3;
cfg.fc_tor_j = [4, 4, 4, 4, 3, 2, 2];
cfg.framelen_sg_j = [0, 0, 0, 0, 21, 21, 21];
fn = fieldnames(opts);
for k = 1:numel(fn), cfg.(fn{k}) = opts.(fn{k}); end
cfg.framelen_sg_j(cfg.framelen_sg_j > 0 & mod(cfg.framelen_sg_j, 2) == 0) = ...
    cfg.framelen_sg_j(cfg.framelen_sg_j > 0 & mod(cfg.framelen_sg_j, 2) == 0) - 1;

[nRow, nJoint] = size(Tor);
fs = 1 / dt;
TorFlt = zeros(nRow, nJoint);

for i = 1:nJoint
    fc = cfg.fc_tor_j(i);
    tor_filt = designfilt('lowpassiir', 'FilterOrder', cfg.filter_order, ...
        'HalfPowerFrequency', fc / (fs / 2), 'DesignMethod', 'butter');

    tor_clean = hampel(Tor(:, i), cfg.hampel_k, cfg.hampel_nsigma);
    tor_lp = filtfilt(tor_filt, tor_clean);

    if cfg.framelen_sg_j(i) > 0
        TorFlt(:, i) = sgolayfilt(tor_lp, cfg.order_sg, cfg.framelen_sg_j(i));
    else
        TorFlt(:, i) = tor_lp;
    end
end

end

%% ========== 对齐裁剪与绘图 ==========

function [idxUse, tCrop, Pos_crop, Vel_crop, VelFlt_crop, Acc_crop, Tor_crop, TorFlt_crop] = ...
        alignAndCropMotionData(Pos, Vel, VelFlt, Acc, Tor, TorFlt, dt, nTrimEdge)

    nRow = size(Pos, 1);
    nAlign = size(Acc, 1);
    idxAlign = 3 : (nRow - 2);
    tAlign = (idxAlign - 1) * dt;

    if 2 * nTrimEdge >= nAlign
        error('trim_edge_s 过大：裁剪后无有效数据，请减小 trim_edge_s');
    end

    idxUse = (1 + nTrimEdge) : (nAlign - nTrimEdge);
    tCrop = tAlign(idxUse);

    Pos_crop = Pos(idxAlign, :);
    Pos_crop = Pos_crop(idxUse, :);
    Vel_crop = Vel(idxAlign, :);
    Vel_crop = Vel_crop(idxUse, :);
    VelFlt_crop = VelFlt(idxAlign, :);
    VelFlt_crop = VelFlt_crop(idxUse, :);
    Acc_crop = Acc(idxUse, :);
    Tor_crop = Tor(idxAlign, :);
    Tor_crop = Tor_crop(idxUse, :);
    TorFlt_crop = TorFlt(idxAlign, :);
    TorFlt_crop = TorFlt_crop(idxUse, :);
end

function plotMotionData(t, Pos, Vel, VelFlt, Acc, Tor, TorFlt)
    colors = lines(7);

    % figure('Name', 'Joint Positon');
    % tiledlayout(7, 1);
    % for i = 1:7
    %     p = nexttile;
    %     plot(p, t, Pos(:, i), 'Color', colors(i, :), 'LineWidth', 1.5);
    %     legend(['P', num2str(i)], 'Location', 'best');
    %     xlabel('t (s)');
    %     ylabel('degree');
    % end
    % 
    % figure('Name', 'Joint Velocity');
    % tiledlayout(7, 1);
    % for i = 1:7
    %     v = nexttile;
    %     plot(v, t, Vel(:, i), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    %     hold on;
    %     plot(v, t, VelFlt(:, i), 'Color', colors(i, :), 'LineWidth', 1.5);
    %     legend('Raw Vel', 'Filtered Vel', 'Location', 'best');
    %     xlabel('t (s)');
    %     ylabel('degree/s');
    %     hold off;
    % end
    % 
    % figure('Name', 'Joint Acceleration');
    % tiledlayout(7, 1);
    % for i = 1:7
    %     a = nexttile;
    %     plot(a, t, Acc(:, i), 'Color', colors(i, :), 'LineWidth', 1.5);
    %     legend('Acc Final', 'Location', 'best');
    %     xlabel('t (s)');
    %     ylabel('degree/s2');
    % end
    % 
    % figure('Name', 'Joint Torque');
    % tiledlayout(7, 1);
    % for i = 1:7
    %     tr = nexttile;
    %     plot(tr, t, Tor(:, i), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    %     hold on;
    %     plot(tr, t, TorFlt(:, i), 'Color', colors(i, :), 'LineWidth', 1.5);
    %     legend('Raw Tor', 'Filtered Tor', 'Location', 'best');
    %     xlabel('t (s)');
    %     ylabel('torque');
    %     hold off;
    % end
end
