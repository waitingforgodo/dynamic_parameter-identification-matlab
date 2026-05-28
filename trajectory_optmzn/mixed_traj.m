function [q,qd,q2d] = mixed_traj(t,C,A,B,w,N)
% ---------------------------------------------------------------------
% 混合轨迹 = 傅里叶轨迹 + 5 阶多项式轨迹
% This function computes "mixed trajectory" meaning the trajectory
% that consistes of finite fourier series and fifth order polynomial
% Inputs:
%   t - time instants at which trajectory should be evaluated
%   C - coefficients of the fifth order polynomail (n_dof x 6)
%   A - coeffincients of the sine in finite fourier series (n_dof x N)
%   B - coefficinets of the cosine in finitne fourier series
%   w - fundamental frequency
%   N - number of harmonics
% ---------------------------------------------------------------------

n = size(A, 1);

[qh,qhd,qh2d] = fourier_series_traj(t,zeros(n,1),A,B,w,N);

qp = C(:,1) + C(:,2).*t + C(:,3).*t.^2 + C(:,4).*t.^3 + ...
                                    C(:,5).*t.^4 + C(:,6).*t.^5;
qpd = C(:,2) + 2*C(:,3).*t + 3*C(:,4).*t.^2 + ...
                                4*C(:,5).*t.^3 + 5*C(:,6).*t.^4;
qp2d = 2*C(:,3) + 6*C(:,4).*t + 12*C(:,5).*t.^2 + 20*C(:,6).*t.^3;

q = qh + qp;
qd = qhd + qpd;
q2d = qh2d + qp2d;

% % ===================== 自动绘图 + 图例 =====================
% figure('Name','混合轨迹波形','Position',[100,100,900,700]);
% 
% % 位置 q
% subplot(3,1,1);
% plot(t, rad2deg(q), 'LineWidth',1.2);
% grid on;
% legendStrings = cell(1, n);
% for i = 1:n
%     legendStrings{i} = sprintf('关节 %d', i);
% end
% legend(legendStrings, 'Location','best'); % 图例
% title('位置 q (°)','fontsize',12);
% xlabel('时间 t (s)');
% ylabel('位置 (°)');
% 
% % 速度 qd
% subplot(3,1,2);
% plot(t, rad2deg(qd), 'LineWidth',1.2);
% grid on;
% legend(legendStrings, 'Location','best');
% title('速度 qd (°/s)','fontsize',12);
% xlabel('时间 t (s)');
% ylabel('速度 (°/s)');
% 
% % 加速度 q2d
% subplot(3,1,3);
% plot(t, rad2deg(q2d), 'LineWidth',1.2);
% grid on;
% legend(legendStrings, 'Location','best');
% title('加速度 q2d (°/s²)','fontsize',12);
% xlabel('时间 t (s)');
% ylabel('加速度 (°/s²)');
% 
% % 总图标题
% sgtitle('混合轨迹：位置 + 速度 + 加速度','fontsize',14,'Color','b');
% % ==================================================================
% 
% end