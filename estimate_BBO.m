%%This script tries to estimate the effect of 7.44m of BBO crystal material
% on the shape of a perfect gaussian pulse

% A gaussian pulse in the time domain is given by E(t) = E0 * exp(-t^2*a)
% The Fourier transform of that pulse is E(w) = E0' * exp(-w^2/(4*a))

c         = 299792458;
thickness = 7.44; % in m
lambda    = 0.980:0.0001:1.080;         % in um
lambda    = lambda * 1e-6;
omega     = 2*pi*c./lambda;
lambda0   = 1.030 * 1e-6;
omega0    = 2*pi*c/lambda0;
f  = c./lambda;
f0 = c./lambda0;

[n, n_0, ~, D] = sellmeier('BBO',lambda0,lambda);
n0 = n{1}; % n(l)
n1 = n{2}; % dn/dl
n2 = n{3}; % d^2n/dl^2

omega2 = linspace(min(omega),max(omega),1000);
lambda2 = 2*pi*c./omega2;
n02 = interp1(omega,n0,omega2, 'pchip');
Phi = omega2./c * thickness .* n02;

%% Analytical formulas for the derivatives of Phi
%  based on the formula for n in sellmeier.m and Phi=w/c*d*n

K = 2*pi*c*1e6;
a = 2.3753;
b = 0.01224;
c = -0.01516;
d = 0.01667;
x = omega2;
ana_GDD = (K^2.*((b.*x.^4.*(3.*d.*x.^2+K^2)+3.*c.*(K^2-d.*x.^2).^3).*(a.*x.^2.*(K^2-d.*x.^2)+b.*x.^4+c.*(K^4-d.*K^2.*x.^2)))./(K^2-d.*x.^2).^4 - K^2.*(c-b./(d-K^2./x.^2).^2).^2) ...
      ./...
      (x.^6.*(a+b./(K^2./x.^2-d)+c*K^2./x.^2).^(3/2));
ana_TOD = -(3.*K^2.*(2.*c.*x.^4.*(K^2-d.*x.^2).*(2.*a^2.*(K^2-d.*x.^2).^5)+a.*b.*x.^2.*(4.*d^4.*x.^8-15.*d^3.*K^2.*x.^6+31.*d^2.*K^4.*x.^4-25.*d.*K^6.*x.^2+5.*K^8)+b^2.*x.^10.*(4.*a^2.*d.*(K^2-d.*x.^2).^2.*(d.*x.^2+K^2)-a.*b.*(8.*d^3.*x.^6-3.*d^2.*K^2.*x.^4-6.*d.*K^4.*x.^2+K^6)+b^2.*d.*x.^4.*(4.*d.*x.^2+K^2))+c^2.*K^2.*x.^2.*(K^2-d.*x.^2).^2.*(5.*a.*(K^2-d.*x.^2).^4+b.*x.^2.*(-5.*d^3.*x.^6+14.*d^2.*K^2.*x.^4-21.*d.*K^4.*x.^2+4.*K^6))+2.*c^3.*K^4.*(K^2-d.*x.^2).^6)) ...
      ./...
      (x.^5.*(K^2-d.*x.^2).^4.*sqrt(a+b./(K^2./x.^2-d)+c.*K^2./x.^2).*(a.*x.^2.*(K^2-d.*x.^2)+b.*x.^4+c.*(K^4-d.*K^2.*x.^2)).^2);
%%

d = omega2(2)-omega2(1);
GD     = diff(Phi) ./ d;
GDD    = diff(GD)  ./ d;
omega3 = linspace(min(omega2),max(omega2(1:end-2)),25);
GDD2   = interp1(omega2(1:end-2),GDD,omega3,'pchip');
TOD    = diff(GDD2) ./ (omega3(2)-omega3(1));

fprintf('GDD: %2.1e fs^2\n',GDD(5)*1e30)
fprintf('TOD: %2.1e fs^3\n',TOD(5)*1e45)

S_w = exp(-(omega2-omega0).^2./1.8e25) .* exp(1i*Phi);

[t,E] = Speck_Fourier(lambda2,S_w);
[m,peak_idx] = max(abs(E));
t_center = t(peak_idx);
center_idx = find_closest_idx(t,0);
E = circshift(E',center_idx-peak_idx)';

phase_minimum = min(Phi);
GDD_phi0 = interp1(omega2(1:end-2),GDD,omega0,'pchip');
TOD_phi0 = interp1(omega3(1:end-1),TOD,omega0,'pchip');
taylor_Phi = 0.5*GDD_phi0*(omega-omega0).^2 + 1/6*TOD_phi0*(omega-omega0).^3;

%% Figure part
%-------------------------figure(1)----------------------------------------
figure(1)
%-------------------------------------subplot(2,2,1)-----------------------
subplot(2,2,1)
[haxes,hline1,hline2] = plotyy(t,abs(E).^2,t,unwrap(angle(E)));
indx1 = find_closest_idx(abs(E).^2,0.5);
indx2 = find_closest_idx(abs(E(1:indx1-10)).^2,0.5);
delta_t = abs(t(indx2) - t(indx1)) * 1e12; % ps
hold on
plot(haxes(1),[t(indx1) t(indx1) t(indx2) t(indx2)],[0 0.5 0.5 0],'r')
hold off
text(1e-12,0.5,sprintf('delta t = %2.2f ps',delta_t),'color','red')
title('Pulse in time with phase curve')
ylabel(haxes(1),'Intensity') % label left y-axis
ylabel(haxes(2),'Phase') % label right y-axis
xlabel(haxes(2),'Time') % label x-axis
axis(haxes(1),[-0.6e-11, 0.6e-11, 0, 1])
axis(haxes(2),[-0.6e-11, 0.6e-11, 50, 150])
set(haxes(1), 'box', 'off');
%-------------------------------------subplot(2,2,2)-----------------------
subplot(2,2,2)
[haxes,hline1,hline2] = plotyy(omega2,abs(S_w).^2,omega,taylor_Phi);
axis(haxes(1),[1.815e15, 1.841e15, 0, 1])
axis(haxes(2),[1.815e15, 1.841e15, 0, 30])
% set(haxes(2),'xlim',[1.815e15, 1.841e15])
%set(haxes(2), 'box', 'off');
title('Spectrum plus the unwrapped phase after the BBO')
ylabel(haxes(1),'Spectrum') % label left y-axis
ylabel(haxes(2),'Phase') % label right y-axis
xlabel(haxes(2),'Angular Frequency') % label x-axis
%-------------------------------------subplot(2,2,3)-----------------------
subplot(2,2,3)
plot(omega2(1:end-2),GDD.*1e30)
hold on
idx0 = find_closest_idx(omega2,omega0);
plot([omega2(idx0) omega2(idx0)],[2.8e5 GDD(idx0)*1e30],'color','red')
plot([1.7e15 omega2(idx0)],[GDD(idx0)*1e30 GDD(idx0)*1e30],'color','red')
hold off
text(1.92e15,3.1e5,sprintf('According to:\nrefractiveindex.info\nthe GVD of BBO at 1.03um is\n42.086 fs^2/mm\n=> 7.44*42.086e3=%2.2fe5fs^2',7.44*0.42086),'EdgeColor','red','FontSize',8,'HorizontalAlignment','right')
text(1.72e15,GDD(idx0)*1e30+0.04e5,sprintf('omega_0 = %2.2e s^{-1}',omega0),'color','red','FontSize',8)
title('Calculation of the GDD')
xlabel('Frequency')
ylabel('GDD in fs^2')
%-------------------------------------subplot(2,2,4)-----------------------
subplot(2,2,4)
plot(omega3(1:end-1),TOD.*1e45)
axis([1.74e15 1.925e15 4.4e5 5.2e5])
title('Calculation of the TOD')
xlabel('Frequency')
ylabel('TOD in fs^3')