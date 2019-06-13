%% SE(3) Visualizer
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                 [Size]
%  T           frame SE(3) to draw                           4*4
%  frame       (optional) existing frame(handler) to update  struct created by this visualizer

%% Outputs
% [Name]      [Description]                                 [Size]
%  frame       frame handler for additional update           struct

%% Examples
% end_effector = plot_SE2(end_effector_T);
% plot_SE2(end_effector_T, end_effector);

%% Implementation
function frame = plot_SE2(T, varargin)
    
    p=T(1:2,3);
    ax=p+T(1:2,1)*0.1;
    ay=p+T(1:2,2)*0.1;
    if (nargin == 1)
        hold on;
        frame.x = plot([p(1),ax(1)], [p(2),ax(2)],'r');
        frame.y = plot([p(1),ay(1)], [p(2),ay(2)],'g');
    elseif(nargin == 2)
        frame = varargin{1};
        set(frame.x, 'XData',[p(1),ax(1)], 'YData',[p(2),ax(2)]);
        set(frame.y, 'XData',[p(1),ay(1)], 'YData',[p(2),ay(2)]);
    end

end