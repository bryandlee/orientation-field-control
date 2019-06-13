%% Analytic logarithm of SO(3)
function w = log_SO3(R)
    if R == eye(3)
        w = zeros(3,1);
    else if 
            
            
    wnorm_sq = w(1)^2 + w(2)^2 + w(3)^2;
    wnorm = wnorm_sq^0.5;
    cw = cos(wnorm);
    sw = sin(wnorm);
    R = [ cw - (w(1)^2*(cw - 1))/wnorm_sq, - (w(3)*sw)/wnorm - (w(1)*w(2)*(cw - 1))/wnorm_sq,   (w(2)*sw)/wnorm - (w(1)*w(3)*(cw - 1))/wnorm_sq;
        (w(3)*sw)/wnorm - (w(1)*w(2)*(cw - 1))/wnorm_sq,                                                   cw - (w(2)^2*(cw - 1))/wnorm_sq, - (w(1)*sw)/wnorm - (w(2)*w(3)*(cw - 1))/wnorm_sq;
        - (w(2)*sw)/wnorm - (w(1)*w(3)*(cw - 1))/wnorm_sq,   (w(1)*sw)/wnorm - (w(2)*w(3)*(cw - 1))/wnorm_sq,                                                   cw - (w(3)^2*(cw - 1))/wnorm_sq];
end


