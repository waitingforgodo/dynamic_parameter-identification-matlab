%%
clc;
clear;
robot = importrobot('Marvin.urdf');

%% 牛顿-欧拉法：
% GenRegNewtonEuler(robot);
%% 拉格朗日法：
% GenRegLagrangian(robot);
% %兼容非正装 
GenRegNewtonEulerGravity(robot);
disp('Done');











































% fileID = fopen('D:\formu.txt', 'w');
% strVar = char(dc{i+1});  
% fprintf(fileID, '%s', strVar);
% fclose(fileID); 




