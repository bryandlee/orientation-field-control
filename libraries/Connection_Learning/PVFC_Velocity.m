function [W_d,V_d,W_d_dot,V_d_dot] = PVFC_Velocity(x,x_dot,R,R_dot,Theta,index,v_bar,v_bar_dot,w_n,delta,eta)
%PVFC_VELOCITY 이 함수의 요약 설명 위치
%   자세한 설명 위치
%% V_d cal
R_g = OrientationFieldGen3D(x,Theta,index);
KK = K(R_g,R,delta);
V_d = v_d(KK,v_bar);
%% W_d cal
k = k_(R_g,R,w_n,eta);
dR_g_v_d = dR_g(x,V_d,Theta,index);
W_d = w_d(R_g,R,k,dR_g_v_d);
%% V_d_dot cal
R_gdot = dR_g(x,x_dot,Theta,index);
DlogR_gR_inv_fro_norm = dlogR_gR_inv_fro_norm(R_g,R,R_gdot,R_dot);
KK_dot = K_dot(KK,DlogR_gR_inv_fro_norm);
V_d_dot = v_d_dot(KK,v_bar,KK_dot,v_bar_dot);
%% W_d_dot cal
k_dot = k_dot_(R_g,R,w_n,DlogR_gR_inv_fro_norm);
W_d_dot = w_d_dot(x,x_dot,R,R_dot,Theta,index,V_d,V_d_dot,k, k_dot,R_g);
end

