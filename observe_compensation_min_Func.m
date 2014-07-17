function val=observe_compensation_min_Func(poly,omega,intensity,phase,lambda,fourierlimit)
    %poly(end-1:end) = 0;
    P = polyval([poly 0 0],omega);
    D = phase - P;
    Sk_cplx = sqrt(intensity) .* exp(1i*D);
    [t,E] = Speck_Fourier(lambda',Sk_cplx);
    t = t * 1e15;
    Int = abs(trapz(t,abs(E).^2));
    E = E .* sqrt(fourierlimit ./ Int);
    val = 1-max(abs(E).^2);
end
% Older version that gave an interesting result
%     P = polyval(x,omega);
%     D = phase - P;
%     Sk_cplx = sqrt(intensity) .* exp(1i*D);
%     [t,E] = Speck_Fourier(lambda'*1e-9,Sk_cplx);
%     t = t * 1e15;
%     Int = abs(trapz(t,abs(E).^2));
%     E = E .* sqrt(fourierlimit ./ Int);
%     val = 1-max(abs(E).^2);

% Good candidate
%     P = polyval(x,omega);
%     D = phase - P;
%     D_start = phase - polyval(start,omega);
%     Sk_cplx = sqrt(intensity) .* exp(1i*D);
%     [t,E] = Speck_Fourier(lambda'*1e-9,Sk_cplx);
%     t = t * 1e15;
%     Int = abs(trapz(t,abs(E).^2));
%     E = E .* sqrt(fourierlimit ./ Int);
%     disp(sum(D_start.^2) / sum(D.^2))
%     val = (1-max(abs(E).^2)) + (sum(D_start.^2) / sum(D.^2))^2 + (sum(D.^2) / sum(D_start.^2))^2;

% Also nice
% val = (1-max(abs(E).^2)) * ((sum(D_start.^2) / sum(D.^2))^2 + (sum(D.^2) / sum(D_start.^2))^2);