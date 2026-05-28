function [T01,T12,T23,T34,T45,T56,T67] = RobotFk(robot,theta1,theta2,theta3,theta4,theta5,theta6,theta7)

    %% 基坐标系 <- J1  
    matB0 = robot.Bodies{1,1}.Joint.JointToParentTransform;
    matT1 = CoordsTrans('R','z',theta1);   
        
    %% J1 <- J2    
    matB1 = robot.Bodies{1,2}.Joint.JointToParentTransform;
    matT2 = CoordsTrans('R','z',theta2);

    %% J2 <- J3  
    matB2 = robot.Bodies{1,3}.Joint.JointToParentTransform;
    matT3 = CoordsTrans('R','z',theta3);

    %% J3 <- J4      
    matB3 = robot.Bodies{1,4}.Joint.JointToParentTransform;
    matT4 = CoordsTrans('R','z',theta4);

    %% J4 <- J5   
    matB4 = robot.Bodies{1,5}.Joint.JointToParentTransform;
    matT5 = CoordsTrans('R','z',theta5);

    %% J5 <- J6   
    matB5 = robot.Bodies{1,6}.Joint.JointToParentTransform;
    matT6 = CoordsTrans('R','z',theta6);

    %% J6 <- J7   
    matB6 = robot.Bodies{1,7}.Joint.JointToParentTransform;
    matT7 = CoordsTrans('R','z',theta7);

    %% 输出 
    T01 = matB0*matT1;
    T12 = matB1*matT2;
    T23 = matB2*matT3;
    T34 = matB3*matT4;
    T45 = matB4*matT5;
    T56 = matB5*matT6;
    T67 = matB6*matT7;

end

%% 坐标系旋转平移变换
function ret = CoordsTrans(type,axis,val)

    if type == 'R'

        % w = val*pi/180;
        w = val;

        if axis == 'x'
            ret = [ 1.0      0.0       0.0       0.0;
                    0.0      cos(w)    -sin(w)   0.0;
                    0.0      sin(w)    cos(w)    0.0;
                    0.0      0.0       0.0       1.0];        
        elseif axis == 'y'
            ret = [ cos(w)   0.0       sin(w)    0.0;
                    0.0      1.0       0.0       0.0;       
                    -sin(w)  0.0       cos(w)    0.0;
                    0.0      0.0       0.0       1.0];
        elseif  axis == 'z'
            ret = [ cos(w)   -sin(w)   0.0       0.0;       
                    sin(w)   cos(w)    0.0       0.0;
                    0.0      0.0       1.0       0.0;
                    0.0      0.0       0.0       1.0];
        else
            ret = 'axis error';
        end

    elseif type == 'P'

        len = val;

        if axis == 'x'
            ret = [ 1.0      0.0       0.0       len;
                    0.0      1.0       0.0       0.0;
                    0.0      0.0       1.0       0.0;
                    0.0      0.0       0.0       1.0];        
        elseif axis == 'y'
            ret = [ 1.0      0.0       0.0       0.0;
                    0.0      1.0       0.0       len;
                    0.0      0.0       1.0       0.0;
                    0.0      0.0       0.0       1.0]; 
        elseif  axis == 'z'
            ret = [ 1.0      0.0       0.0       0.0;
                    0.0      1.0       0.0       0.0;
                    0.0      0.0       1.0       len;
                    0.0      0.0       0.0       1.0]; 
        else
            ret = 'axis error';
        end

    else
        ret = 'type error';
    end

end

