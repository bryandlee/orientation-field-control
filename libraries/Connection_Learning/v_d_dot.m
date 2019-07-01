function outputArg1 = v_d_dot(K,v_bar,K_dot,v_bar_dot)
%V_D_DOT 이 함수의 요약 설명 위치
%   자세한 설명 위치
outputArg1 = K_dot*v_bar+K*v_bar_dot;
end

