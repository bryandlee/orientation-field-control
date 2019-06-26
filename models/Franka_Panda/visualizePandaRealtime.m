%% Visualizer Franka Emika Panda
% 2019 Bryan Dongik Lee

%% Examples
% visualizer = visualizePandaRealtime(robot,q);
% visualizePandaRealtime(robot,q,visualizer);

%% Inputs
% [Name]       [Description]                    [Size]

%% Implementation
function handle = visualizePandaRealtime(robot, q, varargin)
    %% Init
    n = robot.dof;

    T = zeros(4,4,n);
    for i = 1:n
        T(:,:,i) = solveForwardKinematics(q(1:i), robot.A(:,1:i), robot.M(:,:,1:i));
    end

    if nargin == 2
       %% STL Load
        fv_base = stlread(['link0.STL']);  % base link

        for i = 1:n
            handle.fv_zero{i} = stlread(['link' num2str(i) '.STL']);
        end

        handle.fv = handle.fv_zero;
        for i = 1:n
            handle.fv{i}.vertices = (T(1:3,1:3,i)*handle.fv{i}.vertices' + T(1:3,4,i)*ones(1,size(handle.fv{i}.vertices,1)))';
        end
        
       %% Render
        % The model is rendered with a PATCH graphics object. We also add some dynamic
        % lighting, and adjust the material properties to change the specular
        % highlighting.
        handle.fig = figure('Name','Franka Emika Panda','NumberTitle','off','units','pixels','pos',[100 100 800 800]);
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
            handle.render_part{i} = patch(handle.fv{i},'FaceColor',  color(i,:), ...
                     'EdgeColor',       'none',        ...
                     'FaceLighting',    'gouraud',     ...
                     'AmbientStrength', 0.15);
        end

        % draw reference frame
        plot_SE3(eye(4));

        % draw end-effector
        end_effector_T = T(:,:,7) * robot.M_ee;
        handle.end_effector = plot_SE3(end_effector_T);

        % Add a camera light, and tone down the specular highlighting
        camlight('headlight');
        material('dull');
    %     alpha(.5);

        % Fix the axes scaling, and set a nice view angle
        view([-135 35]);
        getframe;
        
    elseif nargin == 3
       %%
        handle = varargin{1};
        for i = 1:n
            handle.fv{i}.vertices = (T(1:3,1:3,i)*handle.fv_zero{i}.vertices' + T(1:3,4,i)*ones(1,size(handle.fv_zero{i}.vertices,1)))';
            set(handle.render_part{i}, 'Vertices', handle.fv{i}.vertices);
        end
%         plot_inertiatensor(T, robot.G);

        end_effector_T = T(:,:,7) * robot.M_ee;
        plot_SE3(end_effector_T, handle.end_effector);
        plot3(end_effector_T(1,4),end_effector_T(2,4),end_effector_T(3,4),'.','Color','r');
        
        getframe;
    end
end