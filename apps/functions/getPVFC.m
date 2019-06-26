%% PVFC
% 

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  qdot         joint vleocities                   n*1
%  G            link inertial matrices             6*6*n
%  Vdot_0       (optional1) base acceleration      6*1

%% Outputs
% [Name]       [Description]                      [Size]

%% Implementation
function tau = getPVFC(A,M,q,qdot,v,vdot,G,Vdot_0,E,M_f,gamma)
    %% Initialization
    n     = size(q,1);          % number of joints
    
    V_0    = zeros(6,1);        % base velocity
 
    T    = zeros(4,4,n);        % T_{i,i-1}
    Ad_T = zeros(6,6,n);        % Ad_T_{i,i-1}

    %% A
    diagA  = zeros(6*n,n);
    for i=1:n
        diagA(6*i-5:6*i,i) = A(:,i);
    end
    
    %% G
    diagG = zeros(6*n,6*n);
    for i=1:n
        diagG(6*i-5:6*i,6*i-5:6*i) = G(:,:,i);
    end
    
    %% ad_V
    V    = zeros(6,n);
    ad_V = zeros(6*n,6*n);
    for i = 1:n
        T(:,:,i)    = exp_se3(-A(:,i)*q(i))*M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        if i == 1
            V(:,i) = Ad_T(:,:,i)*V_0  + A(:,i)*qdot(i);
        else
            V(:,i) = Ad_T(:,:,i)*V(:,i-1)  + A(:,i)*qdot(i);
        end
        
        ad_V(6*i-5:6*i,6*i-5:6*i) = small_ad(V(:,i));
    end

    %% Vdot_base
    Vdot_base = zeros(6*n,1);
    Vdot_base(1:6,1) = Ad_T(:,:,1)* Vdot_0;
    
    %% ad_A_qdot
    ad_A_qdot = zeros(6*n,6*n);
    for i=1:n
        ad_A_qdot(6*i-5:6*i,6*i-5:6*i) = small_ad(A(:,i)*qdot(i));
    end
    
    %% W
    W = zeros(6*n,6*n);
    for i = 2:n
        W(6*i-5:6*i,6*i-11:6*i-6) = Ad_T(:,:,i);
    end
    
    %% L
    L = eye(6*n);
    for i=2:n
        for j=i-1:-1:1
            L(6*i-5:6*i,6*j-5:6*j) = Ad_T(:,:,i)*L(6*i-11:6*i-6,6*j-5:6*j);
        end
    end
    
    %% Dynamic matrices
    M_q = diagA'*L'*diagG*L*diagA;
    c_q = -diagA'*L'*(diagG*L*ad_A_qdot*W + ad_V'*diagG)*L*diagA;
    g_q = diagA'*L'*diagG*L*Vdot_base;
    
    %% Control law
    qdot_bar = zeros(n+1,1);
    qdot_bar(1:n) = qdot;
    qdot_bar(n+1) = sqrt(2*(E - (1/2)*qdot'*M_q*qdot)/M_f);
    
    v_bar = zeros(n+1,1);
    v_bar(1:n) = v;
    v_bar(n+1) = sqrt(2*(E - (1/2)*v'*M_q*v)/M_f);

    vdot_bar = zeros(n+1,1);
    vdot_bar(1:n) = vdot;
    
    % M_dot
    M_q_dot   = zeros(n,n);

    Vdot  = zeros(6,n);
    dVdot = zeros(6,n);
    F     = zeros(6,n);
    dF    = zeros(6,n);
    dtau  = zeros(n,1);

    Ad_T(:,:,n+1) = zeros(6,6);
    F(:,n+1) = zeros(6,1);
    dF(:,n+1) = zeros(6,1);

    for cols = 1:n
        qddot = zeros(n,1);
        qddot(cols) = 1;

        % Vdot
        Vdot(:,1) = A(:,1)*qddot(1);    
        dVdot(:,1) = zeros(6,1);
        for i = 2:n
            T(:,:,i)    = exp_se3(-A(:,i)*q(i))*M(:,:,i);
            Ad_T(:,:,i) = large_Ad(T(:,:,i));
            Vdot(:,i) = Ad_T(:,:,i)*Vdot(:,i-1) + A(:,i)*qddot(i);
            dVdot(:,i) = Ad_T(:,:,i)*dVdot(:,i-1) - small_ad(A(:,i))*Ad_T(:,:,i)*Vdot(:,i-1)*qdot(i);
        end

        % F
        i = n;
        F(:,i) = G(:,:,i)*Vdot(:,i);       
        dF(:,i) = G(:,:,i)*dVdot(:,i);
        for i = n-1:-1:1
            F(:,i) = Ad_T(:,:,i+1)'*F(:,i+1) + G(:,:,i)*Vdot(:,i);
            dF(:,i) = Ad_T(:,:,i+1)'* dF(:,i+1) - Ad_T(:,:,i+1)'*small_ad(A(:,i+1))'*F(:,i+1)*qdot(i+1) + G(:,:,i)*dVdot(:,i);

            dtau(i) = A(:,i)'*dF(:,i);
        end
        M_q_dot(:,cols) = dtau;
    end
    
    vdot_bar(n+1) = -(v'*M_q*vdot + (1/2)*v'*M_q_dot*v)/(v_bar(n+1)*M_f);
    
    
    M_bar = zeros(n+1,n+1);
    M_bar(1:n,1:n) = M_q;
    M_bar(n+1,n+1) = M_f;

    c_bar = zeros(n+1,n+1);
    c_bar(1:n,1:n) = c_q;
    
    g_bar = zeros(n+1,1);
    g_bar(1:n) = g_q;
    
    p = M_bar * qdot_bar;
    P = M_bar * v_bar;
    w = M_bar * vdot_bar + c_bar * v_bar;
    
    tau_c = (w*P' - P*w')*qdot/(2*E);
    tau_f = (P*p' - p*P')*qdot*gamma;

    
    tau = tau_c + tau_f + g_bar;
end