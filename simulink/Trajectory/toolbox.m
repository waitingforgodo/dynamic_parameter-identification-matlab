clc;
clear;
 
%定义连杆的D-H参数
%连杆偏移
d1 = 0;
d2 = 120;
d3 = 0;
d4 = 400;
d5 = 0;
d6 = 0;
%连杆长度
a1 = 0;
a2 = 0;
a3 = 400;
a4 = 10;
a5 = 0;
a6 = 0;
%连杆扭角
alpha1 = 0;
alpha2 = -pi/2;
alpha3 = 0;
alpha4 = -pi/2;
alpha5 = pi/2;
alpha6 = -pi/2;
%建立机器人模型
%       theta  d        a        alpha     
L1=Link([0     d1       a1       alpha1     ],'modified');
L2=Link([0     d2       a2       alpha2     ],'modified');
L3=Link([0     d3       a3       alpha3     ],'modified');
L4=Link([0     d4       a4       alpha4     ],'modified');
L5=Link([0     d5       a5       alpha5     ],'modified');
L6=Link([0     d6       a6       alpha6     ],'modified');
%限制机器人的关节空间
L1.qlim = [(-165/180)*pi,(165/180)*pi];
L2.qlim = [(-150/180)*pi, (60/180)*pi];
L3.qlim = [(-150/180)*pi, (90/180)*pi];
L4.qlim = [(-180/180)*pi,(180/180)*pi];
L5.qlim = [(-115/180)*pi,(115/180)*pi];
L6.qlim = [(-360/180)*pi,(360/180)*pi];
%连接连杆，机器人取名为myrobot
robot=SerialLink([L1 L2 L3 L4 L5 L6],'name','6Rrobot');
robot.base = transl(0 ,0 ,660.4);   %基坐标系进行平移
robot.display();                    %打印出机器人D-H参数表
robot.teach();                      %展示机器人模型
 
%调整坐标轴长度
ax = gca;  
ax.XLim = [-1000, 1000];  
ax.YLim = [-1000, 1000];  
ax.ZLim = [-1000,  1500];  
 
hold on;