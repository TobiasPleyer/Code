function [D2,D3]=PrismS(glas,d,xo)
%     This function calculates the dispersion coefficients for a prism setup
%     like shown below. The origin of the coordinate system is the tip of the
%     second (flipped) prism.
%     Argument explanation:
%             glas: The type of glas/material the prism are made of
%             d   : Distance between first exit and second entry point
%             xo  : Optical path the beam takes through the prisms
%     Variable explanation:
%             lambda  : Wave length in um
%             x_geo   : The geometrical length of the beam in the prisms
% 
%            /\ Prism 1
%           /  \           ________
%          /    \          \      / Prism 2
%         /______\          \    /
%                            \  /
%                             \/
% 
    global A1 A2 A3 l1 l2 l3 w1 w2 w3
    global n
    %---------- SELLMEIER COEFFICIENTS ------------------------------------
    switch glas
        case 'SF11'
            A1 = 1.73759695;
            A2 = 0.313747346;
            A3 = 1.89878101;
            l1 = 0.013188707;
            l2 = 0.0623068142;
            l3 = 155.23629;
            w1 = 2*pi*3/l1*1e14;
            w2 = 2*pi*3/l2*1e14;
            w3 = 2*pi*3/l3*1e14;
    end
    lambda = 1.06; % give lambda in um
    [n,n1,n2,n3] = get_dndl(lambda);
    x_geo = xo/n;
    x = d;
    c = 3e8;

    D2 = lambda^3*1e-6/(pi*c^2) * (-4*x*n1^2 + x_geo*n2 + x_geo*2/lambda*n1) / 1e-30;
    D3 = lambda^4*1e-12/(2*pi^2*c^3) * (12*x * (n1^2 * (1-lambda*n1*(n^-3-2*n)) + ...
         lambda * n1*n2) - x_geo*(3*n2+lambda*n3)) / 1e-45;
end

         function [n,n1,n2,n3]=get_dndl(l)
         global A1 A2 A3 l1 l2 l3
            n = sqrt(1 + (A1 * l^2) / (l^2 - l1) + ...
                         (A2 * l^2) / (l^2 - l2) + ...
                         (A3 * l^2) / (l^2 - l3));
            n1 = (l * (-1) * (A1*l1/(l^2-l1)^2+A2*l2/(l^2-l2)^2+A3*l3/(l^2-l3)^2)) / n;
            n2 = (-10*A1*l^2/(l^2-l1)^2 + 2*A1/(l^2-l1) + 8*A1*l^4/(l^2-l1)^3 + ...
                  -10*A2*l^2/(l^2-l2)^2 + 2*A2/(l^2-l2) + 8*A2*l^4/(l^2-l2)^3 + ...
                  -10*A3*l^2/(l^2-l3)^2 + 2*A3/(l^2-l3) + 8*A3*l^4/(l^2-l3)^3)/(2*n) - ...
                  (2*A1*l/(l^2-l1) - 2*A1*l^3/(l^2-l1)^2 + ...
                   2*A2*l/(l^2-l2) - 2*A2*l^3/(l^2-l2)^2 + ...
                   2*A3*l/(l^2-l3) - 2*A3*l^3/(l^2-l3)^2)/(4*n^3);
            n3 = 3*(2*A1*l/(l^2-l1) - 2*A1*l^3/(l^2-l1)^2 + ...
                    2*A2*l/(l^2-l2) - 2*A2*l^3/(l^2-l2)^2 + ...
                    2*A3*l/(l^2-l3) - 2*A3*l^3/(l^2-l3)^2)^3/(8*n^5) + ...
                    (-24*A1*l/(l^2-l1)^2 - 48*A1*l^5/(l^2-l1)^4 + 72*A1*l^3/(l^2-l1)^3 - ...
                      24*A2*l/(l^2-l2)^2 - 48*A2*l^5/(l^2-l2)^4 + 72*A2*l^3/(l^2-l2)^3 - ...
                      24*A3*l/(l^2-l3)^2 - 48*A3*l^5/(l^2-l3)^4 + 72*A3*l^3/(l^2-l3)^3)/(2*n) - ...
                      3*(-10*A1*l^2/(l^2-l1)^2 + 2*A1/(l^2-l1) + 8*A1*l^4/(l^2-l1)^3 - ...
                          10*A2*l^2/(l^2-l2)^2 + 2*A2/(l^2-l2) + 8*A2*l^4/(l^2-l2)^3 - ...
                          10*A3*l^2/(l^2-l3)^2 + 2*A3/(l^2-l3) + 8*A3*l^4/(l^2-l3)^3) * ...
                          (2*A1*l/(l^2-l1) - 2*A1*l^3/(l^2-l1)^2 + ...
                           2*A2*l/(l^2-l2) - 2*A2*l^3/(l^2-l2)^2 + ...
                           2*A3*l/(l^2-l3) - 2*A3*l^3/(l^2-l3)^2)/(4*n^3);
         end
