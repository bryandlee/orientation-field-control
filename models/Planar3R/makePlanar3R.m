function robot = makePlanar3R()
       
    % degrees of freedom
    robot.dof = 3;

    % link lengths
    l1 = 0.3;
    l2 = 0.2;
    l3 = 0.1;
    
    % screws A_i, i-th screw described in i-th frame
    robot.A = [0 0 0
               0 0 0
               1 1 1
               0 0 0
               0 0 0
               0 0 0];
     
    % link frames M_{i,i-1}
    robot.M(:,:,1) = [1  0  0  0
                      0  1  0  0
                      0  0  1  0
                      0  0  0  1];

    robot.M(:,:,2) = [1  0  0  -l1
                      0  1  0  0
                      0  0  1  0
                      0  0  0  1];
    
    robot.M(:,:,3) = [1  0  0  -l2
                      0  1  0  0
                      0  0  1  0
                      0  0  0  1];    
                  
    robot.M_ee = [1  0  0  l3
                  0  1  0  0
                  0  0  1  0
                  0  0  0  1];
                  
    robot.M_sb = [1  0  0  l1+l2+l3
                  0  1  0  0
                  0  0  1  0
                  0  0  0  1];
    
    % 
    M_si = eye(4,4);
    for i = 1:robot.dof
        M_si = M_si*inverse_SE3(robot.M(:,:,i));
        robot.A_s(:,i) = large_Ad(M_si)*robot.A(:,i);
        robot.A_b(:,i) = large_Ad(inverse_SE3(robot.M_sb))*robot.A_s(:,i);
    end
                        
    % joint limits
%     robot.q_min = [-2.9671 -2.0944 -2.9671 -2.0944 -2.9671 -2.0944 -3.0543]';
%     robot.q_max = [2.9671 2.0944 2.9671 2.0944 2.9671 2.0944 3.0543]';
%     robot.qdot_min = [-1.4835 -1.4835 -1.7453 -1.3090 -2.2689 -2.3562 -2.3562]';
%     robot.qdot_max = [1.4835 1.4835 1.7453 1.3090 2.2689 2.3562 2.3562]';
    
    % inertia matrix Phi for linear dynamics tau = Y * Phi
    for i = 1:robot.dof
        robot.G(:,:,i) = eye(6,6);
        robot.Phi(10*(i-1)+1:10*i) = convertInertiaGToPhi(robot.G(:,:,i));
    end
end