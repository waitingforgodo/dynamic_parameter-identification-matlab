%%
clc;
clear;

addpath(genpath('D:\myWork')); 
robot = importrobot('Marvin.urdf');

%% 牛顿-欧拉法：
[pi_base_NE, baseQR_NE] = GenBaseParameters(robot,'Newton-Euler');

%% 拉格朗日法：
[pi_base_LG, baseQR_LG] = GenBaseParameters(robot,'Lagrangian');

disp('Done');










































% fileID = fopen('D:\formu.txt', 'w');
% strVar = char(dc{i+1});  
% fprintf(fileID, '%s', strVar);
% fclose(fileID); 




