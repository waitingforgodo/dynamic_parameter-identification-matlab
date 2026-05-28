function GenRegNewtonEulerGravity(robot)
%GENREGNEWTONEULERGRAVITY 含可变重力方向的 Newton-Euler 回归矩阵
%
%   GenRegNewtonEulerGravity(robot)
%
%   观测方程:  tau = Y(q, dq, ddq, g) * pi
%   g          基座坐标系下的重力线加速度 [3×1] (m/s^2)
%              标准竖装时约为 [0; 0; 9.8065]，侧装/倒装由基座姿态换算
%
%   输出: autogen/standard_regressor_marvin.m
%         Y = standard_regressor_marvin(q, dq, ddq, g)

    %% 符号
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nJoints = 7;
    q_sym = sym('q%d', [nJoints, 1], 'real');
    dq_sym = sym('dq%d', [nJoints, 1], 'real');
    ddq_sym = sym('ddq%d', [nJoints, 1], 'real');
    g_sym = sym('g', [3, 1], 'real');   % 基座系重力线加速度 g

    w0 = sym(zeros(3, 1));
    dw0 = sym(zeros(3, 1));
    dv0 = g_sym;

    %% 计算  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 运动学Fk  
    [T01,T12,T23,T34,T45,T56,T67] = RobotFk(robot,q_sym(1),q_sym(2),q_sym(3),q_sym(4),q_sym(5),q_sym(6),q_sym(7));  

    R01 = T01(1:3,1:3); P01 = T01(1:3,4);
    R12 = T12(1:3,1:3); P12 = T12(1:3,4);
    R23 = T23(1:3,1:3); P23 = T23(1:3,4);
    R34 = T34(1:3,1:3); P34 = T34(1:3,4);
    R45 = T45(1:3,1:3); P45 = T45(1:3,4);
    R56 = T56(1:3,1:3); P56 = T56(1:3,4);
    R67 = T67(1:3,1:3); P67 = T67(1:3,4);

    % 牛顿-欧拉 外推公式 
    [w1,dw1,dv1] = GenerateLinkSpeed('R',w0,dw0,dv0,R01,P01,dq_sym(1),ddq_sym(1)); 
    [w2,dw2,dv2] = GenerateLinkSpeed('R',w1,dw1,dv1,R12,P12,dq_sym(2),ddq_sym(2)); 
    [w3,dw3,dv3] = GenerateLinkSpeed('R',w2,dw2,dv2,R23,P23,dq_sym(3),ddq_sym(3)); 
    [w4,dw4,dv4] = GenerateLinkSpeed('R',w3,dw3,dv3,R34,P34,dq_sym(4),ddq_sym(4)); 
    [w5,dw5,dv5] = GenerateLinkSpeed('L',w4,dw4,dv4,R45,P45,dq_sym(5),ddq_sym(5)); 
    [w6,dw6,dv6] = GenerateLinkSpeed('L',w5,dw5,dv5,R56,P56,dq_sym(6),ddq_sym(6));
    [w7,dw7,dv7] = GenerateLinkSpeed('L',w6,dw6,dv6,R67,P67,dq_sym(7),ddq_sym(7));

    [A1] = BuildMatA(dv1,dw1,w1);
    [A2] = BuildMatA(dv2,dw2,w2);
    [A3] = BuildMatA(dv3,dw3,w3);
    [A4] = BuildMatA(dv4,dw4,w4);
    [A5] = BuildMatA(dv5,dw5,w5);
    [A6] = BuildMatA(dv6,dw6,w6);
    [A7] = BuildMatA(dv7,dw7,w7);

    T12 = [R12 zeros(3,3); STf(P12)*R12 R12];
    T23 = [R23 zeros(3,3); STf(P23)*R23 R23];
    T34 = [R34 zeros(3,3); STf(P34)*R34 R34];
    T45 = [R45 zeros(3,3); STf(P45)*R45 R45];
    T56 = [R56 zeros(3,3); STf(P56)*R56 R56];
    T67 = [R67 zeros(3,3); STf(P67)*R67 R67];

    U11 = A1; U12 = T12*A2; U13 = T12*T23*A3; U14 = T12*T23*T34*A4; U15 = T12*T23*T34*T45*A5; U16 = T12*T23*T34*T45*T56*A6; U17 = T12*T23*T34*T45*T56*T67*A7;
    U22 = A2; U23 = T23*A3; U24 = T23*T34*A4; U25 = T23*T34*T45*A5; U26 = T23*T34*T45*T56*A6; U27 = T23*T34*T45*T56*T67*A7;
    U33 = A3; U34 = T34*A4; U35 = T34*T45*A5; U36 = T34*T45*T56*A6; U37 = T34*T45*T56*T67*A7;
    U44 = A4; U45 = T45*A5; U46 = T45*T56*A6; U47 = T45*T56*T67*A7;
    U55 = A5; U56 = T56*A6; U57 = T56*T67*A7;
    U66 = A6; U67 = T67*A7;
    U77 = A7;

    U = [       U11           U12           U13           U14           U15           U16           U17;
         zeros(6,10)          U22           U23           U24           U25           U26           U27;
         zeros(6,10)   zeros(6,10)          U33           U34           U35           U36           U37;
         zeros(6,10)   zeros(6,10)   zeros(6,10)          U44           U45           U46           U47;
         zeros(6,10)   zeros(6,10)   zeros(6,10)   zeros(6,10)          U55           U56           U57;
         zeros(6,10)   zeros(6,10)   zeros(6,10)   zeros(6,10)   zeros(6,10)          U66           U67;
         zeros(6,10)   zeros(6,10)   zeros(6,10)   zeros(6,10)   zeros(6,10)   zeros(6,10)          U77];

    Y = [U(6,:); 
         U(12,:); 
         U(18,:);
         U(24,:); 
         U(30,:);
         U(36,:); 
         U(42,:)];

    % Yi = simplify(Y);

    % 生成可高速调用的回归矩阵函数（与 GenRegNewtonEuler 一致，输出到 autogen/）
    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    autogenDir = fullfile(repoRoot, 'autogen');
    if ~exist(autogenDir, 'dir') && ~mkdir(autogenDir)
        error('GenRegNewtonEulerGravity:MkdirFailed', '无法创建 autogen 目录: %s', autogenDir);
    end
    outFile = fullfile(autogenDir, 'standard_regressor_marvin');
    matlabFunction(Y, 'File', outFile, 'Vars', {q_sym, dq_sym, ddq_sym, g_sym});

    fprintf('NE 重力回归已生成: %s.m\n', outFile);
    fprintf('  调用: Y = standard_regressor_marvin(q, dq, ddq, g)\n');

end

%% 生成杆件的速度 
function [w,dw,dv] = GenerateLinkSpeed(type,w0,dw0,dv0,R,P,dq,ddq)
% 计算杆件的角速度、角加速度和线加速度 
    Z = [0.0 0.0 1.0]';
    if type == 'R'
        w = R'*w0 + dq*Z;
        dw = R'*dw0 + cross(R'*w0,dq*Z) + ddq*Z;
        dv = R'*(cross(dw0,P) + cross(w0,cross(w0,P)) + dv0);
    elseif type == 'L'
        w = R'*w0;
        dw = R'*dw0;
        dv = R'*(cross(dw0,P) + cross(w0,cross(w0,P)) + dv0) + 2.0*(cross(w,dq*Z)) + ddq*Z;  
    end
end

%% 建立A矩阵 
function [retA] = BuildMatA(dv,dw,w)
    matA(1:3,1) = dv;
    matA(1:3,2:4) = STf(dw) + STf(w)* STf(w);
    matA(1:3,5:10) = 0;
    
    matA(4:6,1) = 0;
    matA(4:6,2:4) = -STf(dv);
    matA(4:6,5:10) = KTf(dw) + STf(w)*KTf(w);
    
    retA = matA;
end

%% S变换
function [ret] = STf(vec)
    ret = [     0    -vec(3)    vec(2); 
            vec(3)        0    -vec(1); 
           -vec(2)    vec(1)        0];
end

%% K变换
function [ret] = KTf(vec)
    ret = [ vec(1)   vec(2)   vec(3)       0        0        0  ; 
                0    vec(1)       0    vec(2)   vec(3)       0  ; 
                0        0    vec(1)       0    vec(2)   vec(3)];
end
