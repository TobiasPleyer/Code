function [D2,D3]=Prism(glas,d,xo)
%     This function calculates the dispersion coefficients for a prism setup
%     like shown below. The origin of the coordinate system is the tip of the
%     second (flipped) prism.
%     Argument explanation:
%             glas: The type of glas/material the prism are made of
%             d   : Distance between first exit and second entry point
%             xo  : Optical path the beam takes through the prisms
%     Variable explanation:
%             lambda  : Wave length in um
%             omega   : The angular frequency equivalent to lambda
%             gammaB  : The Brewster angle for the respective central wavelength
%             alphaA  : The Apex angle of the prism
%             e_prism : 
%             xo      : The optical length of the beam inside the prisms
%             x_geo   : The geometrical length of the beam in the prisms
% 
%            /\ Prism 1
%           /  \           ________
%          /    \          \      / Prism 2
%         /______\          \    /
%                            \  /
%                             \/
% 
    global A1 A2 A3 l1 l2 l3 w1 w2 w3 alphaA lambda omega
    global gammaB n %e_prism x_w0 y_w0
%     glas = 'SF11';
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
    omega = 2*pi*3/lambda*1e14;
    [n,n1,n2,n3] = get_dndl(lambda);
    gammaB = atan(n);
    alphaA = 2*asin(cos(gammaB));
%     e_prism = 0.002;
%     x_w0 = -0.09; % x-coordinate where beam leaves first prism
%     y_w0 = 0.06; % y-coordinate where beam leaves first prism
%     x_w0 = x;
%     y_w0 = y;
%     xo = optical_path(omega);
    x_geo = xo/n;
    x = d;
%     x = sqrt(x_w0^2 + y_w0^2);
    c = 3e8;
%     [e,e1,e2,e3,n,n1,n2]=get_all_derivatives(omega);
%     fprintf('xo: %g\n',xo)
%     fprintf('e: %g\n',e)
%     fprintf('e1: %g\n',e1)
%     fprintf('e2: %g\n',e2)
%     fprintf('e3: %g\n',e3)
%     fprintf('n: %g\n',n)
%     fprintf('n1: %g\n',n1)
%     fprintf('n2: %g\n',n2)
%     fprintf('w1: %g\n',w1)
%     fprintf('w2: %g\n',w2)
%     fprintf('w3: %g\n',w3)
%     fprintf('alphaA: %g\n',alphaA)
%     fprintf('omega: %g\n',omega)
%     fprintf('gammaB: %g\n',gammaB)
%     D2 = -4*lambda*1e-6/(2*pi*c^2)*n1^2*x + xo/c * e2;
%     D3 = 12*x*(lambda*1e-6)^4/(4*pi^2*c^3) * ...
%          (n1^2 * (1 - lambda*1e-6*n1 * (1/n - 2*n)) + ...
%          lambda*1e-6 * (n1 * n2)) + xo/c * e3;
%     fprintf('n: %g\n',n)
%     fprintf('n1: %g\n',n1)
%     fprintf('n2: %g\n',n2)
%     fprintf('n3: %g\n',n3)
    D2 = lambda^3*1e-6/(pi*c^2) * (-4*x*n1^2 + x_geo*n2 + x_geo*2/lambda*n1) / 1e-30;
    D3 = lambda^4*1e-12/(2*pi^2*c^3) * (12*x * (n1^2 * (1-lambda*n1*(n^-3-2*n)) + ...
         lambda * n1*n2) - x_geo*(3*n2+lambda*n3)) / 1e-45;
end
         
%          function n=n_l(lambda)
%          global A1 A2 A3 l1 l2 l3
%             n = sqrt(1 + (A1 * lambda^2) / (lambda^2 - l1) + ...
%                          (A2 * lambda^2) / (lambda^2 - l2) + ...
%                          (A3 * lambda^2) / (lambda^2 - l3));
%          end
%          
%          function n=n_w(omega)
%          global A1 A2 A3 l1 l2 l3 lambda
%             lambda = 2*pi * 3e14 / omega;
%             n = sqrt(1 + (A1 * lambda^2) / (lambda^2 - l1) + ...
%                          (A2 * lambda^2) / (lambda^2 - l2) + ...
%                          (A3 * lambda^2) / (lambda^2 - l3));
%          end
%          % Definition out of Karsch script to calculate the dispersion
%          % coefficients
%          function e=eta_w(omega)
%             e = n_w(omega) * omega;
%          end
%          % First derivative
%          function e=eta1_w(omega)
%          global A1 A2 A3 w1 w2 w3
%             e = 1/n_w(omega) + (1/n_w(omega)) * ((A1 * w1^4) / (w1^2-omega^2)^2 + ...
%                                                  (A2 * w2^4) / (w2^2-omega^2)^2 + ...
%                                                  (A3 * w3^4) / (w3^2-omega^2)^2);
%          end
%          % Second derivative
%          function e=eta2_w(omega)
%          global A1 A2 A3 w1 w2 w3
%             e = eta1_w(omega) * (1/omega - eta1_w(omega)/eta_w(omega)) + ...
%                 4 * omega^2/eta_w(omega) * ((A1 * w1^4) / (w1^2-omega^2)^3 + ...
%                                             (A2 * w2^4) / (w2^2-omega^2)^3 + ...
%                                             (A3 * w3^4) / (w3^2-omega^2)^3);
%          end
%          % Third derivative
%          function e=eta3_w(omega)
%          global A1 A2 A3 w1 w2 w3
%             e = 3*(eta2_w(omega)-eta1_w(omega)/omega) * ...
%                 (1/omega-eta1_w(omega)/eta_w(omega)) + ...
%                 24*omega^3/eta_w(omega) * ((A1 * w1^4) / (w1^2-omega^2)^4 + ...
%                                            (A2 * w2^4) / (w2^2-omega^2)^4 + ...
%                                            (A3 * w3^4) / (w3^2-omega^2)^4);
%          end
%          
%          function [e,e1,e2,e3]=get_eta(omega)
%          global A1 A2 A3 w1 w2 w3
%             e  = n_w(omega) * omega;
%             e1 = omega/e + (omega/e) * ((A1 * w1^4) / (w1^2-omega^2)^2 + ...
%                                         (A2 * w2^4) / (w2^2-omega^2)^2 + ...
%                                         (A3 * w3^4) / (w3^2-omega^2)^2);
%             e2 = e1 * (1/omega - e1/e) + ...
%                  4 * omega^2/e * ((A1 * w1^4) / (w1^2-omega^2)^3 + ...
%                                   (A2 * w2^4) / (w2^2-omega^2)^3 + ...
%                                   (A3 * w3^4) / (w3^2-omega^2)^3);
%             e3 = 3*(e2-e1/omega) * ...
%                  (1/omega-e1/e) + ...
%                  24*omega^3/e * ((A1 * w1^4) / (w1^2-omega^2)^4 + ...
%                                  (A2 * w2^4) / (w2^2-omega^2)^4 + ...
%                                  (A3 * w3^4) / (w3^2-omega^2)^4);
%          end
%          
%          function [e,e1,e2,e3,n,n1,n2]=get_all_derivatives(omega)
%             [e,e1,e2,e3] = get_eta(omega);
%             n = n_w(omega);
%             n1 = (1/omega) * (e1-n);
%             n2 = (1/omega) * (e2-2*n1);
%          end
         
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
         
%          function xo=optical_path(omega)
%          global alphaA e_prism gammaB x_w0 y_w0 n
%             g1 = (2*e_prism) * sin(alphaA/2);
% %             fprintf('g1: %g\n',g1)
%             gamma_w = asin(n*sin(alphaA - ...
%                       asin(sin(gammaB)/n)));
% %             fprintf('gamma_w: %g\n',gamma_w)
%             b = g1;
%             delta2 = alphaA - asin(sin(gammaB)/n);
%             delta1 = alphaA - delta2;
%             %fprintf('delta1: %g\n',delta1)
%             %fprintf('delta2: %g\n',delta2)
%             beta  = 0.5*pi - delta2;
%             alpha = alphaA/2 - delta1;
%             gamma = pi - alpha - beta;
%             c = b * sin(gamma) / sin(beta);
%             a = b * sin(alpha) / sin(beta);
%             %fprintf('a: %g\n',a)
%             x_w = x_w0 + a * sin(alphaA/2);
%             y_w = y_w0 - a * cos(alphaA/2);
% %             fprintf('x_w: %g\n',x_w)
% %             fprintf('y_w: %g\n',y_w)
%             x_eprime = 1/(tan(gamma_w-alphaA/2)- ...
%                           tan(0.5*pi-alphaA/2)) * ...
%                           (tan(gamma_w-alphaA/2)*x_w + y_w);
%             y_eprime = -tan(0.5*pi-alphaA/2) * x_eprime;
%             fprintf('x_eprime: %g\n',x_eprime)
%             fprintf('y_eprime: %g\n',y_eprime)
% %             fprintf('sigma: %g\n',tan(gamma_w-alphaA/2))
% %             fprintf('origin: %g\n',tan(0.5*pi-alphaA/2))
%             e_prime = sqrt(x_eprime^2 + y_eprime^2);
% %             fprintf('eprime: %g\n',e_prime)
%             f = e_prism + a;
%             g2 = e_prime / f * g1;
% %             fprintf('g2: %g\n',g2)
%             x = c + g2;
%             xo = n * x;
%          end
