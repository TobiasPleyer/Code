function [D2,D3]=Grating(D0)
    c = 3e8;
    N = 300/1e-3; % the grating groove density
    lambda = 1.06e-6;
    gamma = 8.6/360*2*pi;
    beta = asin(N*lambda-sin(gamma));
    D = D0/cos(beta);
    D2 = -2*D/c*lambda/(2*pi*c) * (lambda*N/cos(beta))^2;
    D3 = -3*D2*lambda/(2*pi*c) * (1 + lambda*N*sin(beta)/(cos(beta))^2);
    D2 = D2 / 1e-30;
    D3 = D3 / 1e-45;
%     D4 = 3*D2*(lambda/2/pi/c)^2 * ...
%          (4*(1 + lambda*N*sin(beta)/(cos(beta))^2) + ...
%           (lambda*N/(cos(beta))^2)^2);
end