%% Visualizer Franka Emika Panda
% 2019 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                    [Size]
%  trajectory   robot trajectory to visualize    7*num_time

%% Implementation
function visualizePanda(varargin)
    %% Robot Init
    robot = makeFrankaPanda();

    export_video = false;
       
    n = robot.dof;
    q = zeros(n,1);
    T = zeros(4,4,n);
    for i = 1:n
        T(:,:,i) = solveForwardKinematics(q(1:i), robot.A(:,1:i), robot.M(:,:,1:i));
    end

    trajectory = rand(7,1);
    if nargin > 0
        if size(varargin{1},1) == n
            trajectory = varargin{1};
        else
            error('number of the joints and the trajectory does not match');
        end
    end
    
    num_time = size(trajectory, 2);
    
    %% STL Load
    fv_base = stlread(['link0.STL']);  % base link

    for i = 1:n
        fv_zero{i} = stlread(['link' num2str(i) '.STL']);
    end

    fv = fv_zero;
    for i = 1:n
        fv{i}.vertices = (T(1:3,1:3,i)*fv{i}.vertices' + T(1:3,4,i)*ones(1,size(fv{i}.vertices,1)))';
    end

    %% Video
    if export_video
        writerObj = VideoWriter('kuka','MPEG-4'); % 
        writerObj.FrameRate = 30;
    
        % open the video writer
        open(writerObj);
    end
    
    %% Render
    % The model is rendered with a PATCH graphics object. We also add some dynamic
    % lighting, and adjust the material properties to change the specular
    % highlighting.
    figure('Name','Franka Emika Panda','NumberTitle','off','units','pixels','pos',[100 100 800 800]);
    hold on;
    axis equal;
    axis([-1 1 -1 1 -0.5 1.3]);
    xlabel('x'); ylabel('y'); zlabel('z');

    % draw ground
%     ground = [-1,-1, 1, 1;
%               -1, 1, 1,-1;
%                0, 0, 0, 0]*2;
%     fill3(ground(1,:),ground(2,:),ground(3,:),'g'); alpha(0.2);

    % draw base link
    patch(fv_base,'FaceColor',       [0.7 0.7 0.7], ...
                 'EdgeColor',       'none',        ...
                 'FaceLighting',    'gouraud',     ...
                 'AmbientStrength', 0.15);

    % draw 7 links
    color = ones(n,3)*0.8;
    color(2,:) = [100 100 100]/255;
    color(6,:) = [100 100 100]/255;

    for i = 1:7
        render_part{i} = patch(fv{i},'FaceColor',  color(i,:), ...
                 'EdgeColor',       'none',        ...
                 'FaceLighting',    'gouraud',     ...
                 'AmbientStrength', 0.15);
    end
    
    % draw reference frame
    plot_SE3(eye(4));
    
    % draw end-effector
    end_effector_M = robot.M_ee;
    end_effector_T = T(:,:,7) * end_effector_M;
    end_effector = plot_SE3(end_effector_T);
    
    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');
%     alpha(.5);

    % Fix the axes scaling, and set a nice view angle
    view([-135 35]);
    getframe;
    
%     T_d = zeros(4,4,num_time);
    %%
    trajectory_plot = zeros(3,num_time);
    for time = 1:num_time
        for i = 1:n
            T(:,:,i) = solveForwardKinematics(trajectory(1:i,time), robot.A(:,1:i), robot.M(:,:,1:i));
            fv{i}.vertices = (T(1:3,1:3,i)*fv_zero{i}.vertices' + T(1:3,4,i)*ones(1,size(fv_zero{i}.vertices,1)))';
            set(render_part{i}, 'Vertices', fv{i}.vertices);
        end
%         plot_inertiatensor(T, robot.G);

        end_effector_T = T(:,:,7) * end_effector_M;
        trajectory_plot(:,time)=end_effector_T(1:3,4);
        
%         T_d(:,:,time) = end_effector_T;
    end
    
    plot3(trajectory_plot(1,:),trajectory_plot(2,:),trajectory_plot(3,:),'.','Color','r');

    
    %% Animation Loop
    for iter =1:1
        for time = 1:num_time
            for i = 1:n
                T(:,:,i) = solveForwardKinematics(trajectory(1:i,time), robot.A(:,1:i), robot.M(:,:,1:i));
                fv{i}.vertices = (T(1:3,1:3,i)*fv_zero{i}.vertices' + T(1:3,4,i)*ones(1,size(fv_zero{i}.vertices,1)))';
                set(render_part{i}, 'Vertices', fv{i}.vertices);
            end

            end_effector_T = T(:,:,7) * end_effector_M;
            plot_SE3(end_effector_T, end_effector);

            frame = getframe(gcf);
            
            if export_video
                writeVideo(writerObj, frame);
            end
        end
    end
    
    if export_video
        % close the writer object
        close(writerObj);
    end
end