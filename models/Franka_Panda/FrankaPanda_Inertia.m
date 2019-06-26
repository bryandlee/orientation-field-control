phi = zeros(10,7);

phi(1,1) = 2.3599995791;
phi(5,1) = 0.0186863903;
phi(6,1) = 0.0143789874;
phi(7,1) = 0.00906812;

phi(1,2) = 2.379518833;
phi(5,2) = 0.0190388734;
phi(6,2) = 0.0091429124;
phi(7,2) = 0.014697537;

phi(1,3) = 2.6498823337;
phi(5,3) = 0.0129300178;
phi(6,3) = 0.0150242121;
phi(7,3) = 0.0142734598;

phi(1,4) = 2.6948018744;
phi(5,4) = 0.0133874611;
phi(6,4) = 0.014514325;
phi(7,4) = 0.01551755518;

phi(1,5) = 2.9812816864;
phi(5,5) = 0.0325565705;
phi(6,5) = 0.0270660472;
phi(7,5) = 0.0115023375;

phi(1,6) = 1.1285806309;
phi(5,6) = 0.0026052565;
phi(6,6) = 0.0039897229;
phi(7,6) = 0.0047048591;

phi(1,7) = 0.4052912465;
phi(5,7) = 0.0006316592;
phi(6,7) = 0.0006319639;
phi(7,7) = 0.0010607721;  

for i = 1:7

    fv{i} = stlread(['link' num2str(i) '.STL']);
 
    figure(i)
    hold on;
    axis equal;
    axis([-1 1 -1 1 -0.5 1.1]);
    xlabel('x'); ylabel('y'); zlabel('z');

    plot_SE3(eye(4));
    
    patch(fv{i},'FaceColor',       [0.7 0.7 0.7], ...
                 'EdgeColor',       'none',        ...
                 'FaceLighting',    'gouraud',     ...
                 'AmbientStrength', 0.15); alpha(0.5);

    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');

    % Fix the axes scaling, and set a nice view angle
    view([-135 35]);
    getframe;
    
    
    centroid = mean(fv{i}.vertices);
    plot3(centroid(1),centroid(2),centroid(3), '.r');
    
    G(:,:,i) = convertInertiaPhiToG(phi(:,i));
    T = eye(4);
    T(1:3,4) = centroid';
    T = inverse_SE3(T);
    G(:,:,i) = large_Ad(T)'*G(:,:,i)*large_Ad(T);
    
    plot_inertiatensor(eye(4), G(:,:,i));
end

Phi = zeros(70,1);
for i = 1:7
    Phi(10*(i-1)+1:10*i) = convertInertiaGToPhi(G(:,:,i));
end


