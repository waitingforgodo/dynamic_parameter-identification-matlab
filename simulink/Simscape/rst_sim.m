%%
clc;
clear;
close all;
%%
addpath('../URDF/XM7p/urdf/');
addpath('../URDF');
XM7p = importrobot('XM7p.urdf');
showdetails(XM7p);
show(XM7p,'Frames','on','Visuals','on');
XM7p.Gravity = [0 0 -9.81];
%%
config = homeConfiguration(XM7p);
show(XM7p)
config(1).JointPosition = pi/4;
config(2).JointPosition = pi/6;
show(XM7p,config);
%% load urdf in simulink
XM7p_SC = smimport('XM7p.urdf');