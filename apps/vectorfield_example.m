%% Vector Field
x1 = 0:0.05:1;
x2 = x1;
x3 = x1;
[x1, x2, x3] = meshgrid(x1, x2, x3);
x1dot = - (x1 - 0.5);
x2dot = - (x2 - 0.5);
x3dot = - (x3 - 0);

% quiver(x1, x2, x1dot, x2dot);
startx1 = [linspace(0,0,10), linspace(0,1,10), linspace(1,1,10), linspace(1,0,10)];
startx2 = [linspace(0,1,10), linspace(1,1,10), linspace(1,0,10), linspace(0,0,10)];
streamslice(x1,x2,x3,x1dot,x2dot,x3dot,[],[],0);