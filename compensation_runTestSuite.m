clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
format shortg;
warning off
%format compact;
set(0,'DefaultFigureWindowStyle','docked')


fprintf('---------RUNNING TEST SUITE---------\n')
fprintf('l0: 1030 nm\n')
fprintf('w : 1.82:0.0001:1.85 * 1e15 s^-1\n')
fprintf('I : exp(-(w-w0).^2/1e-5)\n')
fprintf('p : 80000*(w-w0).^2\n')
fprintf('\n')
figNum = 1;
c  = 299792458;
l0 = 1030;
w0 = 1e-15*  2*pi*c / (l0*1e-9);

% TEST 0: Starting conditions
fprintf('\n')
fprintf('TEST 0: Starting conditions\n')

w_low  = 1.82;
w_high = 1.84;
w_inc  = 0.0001;
w  = w_low:w_inc:w_high;
l  = 1e9*    2*pi*c ./ (w*1e15);
l  = fliplr(l);
p  = 80000*(w-w0).^2;
I  = exp(-(w-w0).^2/1e-5);

[Int_F,t_F,Ek_F] = compensation_calcFourierlimit(I,l);
[Int,t,E]        = compensation_calcFourierlimit(I,l,p);

figure(figNum)
    figNum = figNum + 1;
    [AX,H1,H2] = plotyy(w,I,w,p);
    set(AX,'xlim',[w_low w_high]);
    set(get(AX(1),'Ylabel'),'String','arbitrary units')
    set(get(AX(2),'Ylabel'),'String','[rad]')
    xlabel('omega [1/fs]')

% TEST 1: See how our Fourier transformation performs.
fprintf('TEST 1: Fourier Trafo\n')

p1  = 50000*(w-w0).^2;
p2  = 80000*(w-w0).^2;
p3  = 110000*(w-w0).^2;
[Int1,t1,E1] = compensation_calcFourierlimit(I,l,p1);
[Int2,t2,E2] = compensation_calcFourierlimit(I,l,p2);
[Int3,t3,E3] = compensation_calcFourierlimit(I,l,p3);
factor = Int_F/Int;

figure(figNum)
    figNum = figNum + 1;
    plot(t_F,abs(Ek_F).^2,'g')
    hold on
    plot(t1,abs(E1).^2.*factor,'r')
    plot(t2,abs(E2).^2.*factor,'b')
    plot(t3,abs(E3).^2.*factor,'m')
    hold off
    xlim([-2000 2000])
    xlabel('time [fs]')
    ylabel('relative units')
    legend('Fourier limit','Real pulse GDD 5e4','Real pulse GDD 8e4','Real pulse GDD 11e4')
    title('Observation of the effects of different GDDs on the pulse')
    
% TEST 2: See if adding zeros increases our sample rate.
fprintf('TEST 2: Increase sample rate\n')

N = length(w);
[w_ext,I_ext,p_ext,l_ext] = compensation_extendPhaseByZeros(w,I,p,N);
l_ext = fliplr(l_ext);
[Int_F2,t_F2,Ek_F2] = compensation_calcFourierlimit(I_ext,l_ext);
[Int2,t2,E2] = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
factor2 = Int_F2/Int2;

N = 2*length(w);
[w_ext,I_ext,p_ext,l_ext] = compensation_extendPhaseByZeros(w,I,p,N);
l_ext = fliplr(l_ext);
[Int_F3,t_F3,Ek_F3] = compensation_calcFourierlimit(I_ext,l_ext);
[Int3,t3,E3] = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
factor3 = Int_F3/Int3;

figure(figNum)
    figNum = figNum + 1;
    plot(t_F,abs(Ek_F).^2,'g')
    hold on
    plot(t,abs(E).^2.*factor,'r')
    plot(t_F2,abs(Ek_F2).^2,'b')
    plot(t2,abs(E2).^2.*factor2,'c')
    plot(t_F3,abs(Ek_F3).^2,'m')
    plot(t3,abs(E3).^2.*factor3,'k')
    hold off
    xlim([-2000 2000])
    xlabel('time [fs]')
    ylabel('relative units')
    legend('Fourier limit 1xL','Real pulse 1xL','Fourier limit 2xL','Real pulse 2xL','Fourier limit 3xL','Real pulse 3xL')
    title('Observation of the increase in peak sampling rate by adding zeros')
    
% TEST 3: See if our brute force ansatz works with our simple sample phase
fprintf('TEST 3: Approximation of the phase curve via brute force\n')

fprintf('------> 1D\n')
p = 80000*(w-w0).^2;
% Reserve arrays for storage
min_a = 10000;
max_a = 100000;
n_a   = 10;
A = linspace(min_a,max_a,n_a);
V = zeros(1,n_a);

% Enter a loop to brute force the solutions
N = 1;
% CAREFUL: OUR POLYNOMIAL IS 1e4(w^2 - 2ww0 + w0^2) --> [8e4 -3.6576*8e4 3.3445*8e4]
b = -2*w0*8e4;
c = w0^2*8e4;
figure(5)
hold on
for a=A
    poly = [a b c];
    plot(w,polyval(poly,w))
    % Calculate the value
    val = compensation_minFuncForBruteForce_no_global(poly,w,I,p,Int_F);
    V(N) = val;
    if mod(N,1000)==0
        fprintf('N: %d\n',N)
    end
    N = N + 1;
end
hold off

figure(figNum)
    figNum = figNum + 1;
    plot(V)

fprintf('------> 2D\n')
% p  = 5000*(w-w0).^3+80000*(w-w0).^2;
% figure(figNum)
%     figNum = figNum + 1;
%     [AX,H1,H2] = plotyy(w,I,w,p);
%     set(AX,'xlim',[w_low w_high]);
%     set(get(AX(1),'Ylabel'),'String','arbitrary units')
%     set(get(AX(2),'Ylabel'),'String','[rad]')
%     xlabel('omega [1/fs]')
% 
% % Reserve arrays for storage
% A = linspace(min_a,max_a,n_a);
% B = linspace(min_b,max_b,n_b);
% V = zeros(n_b,n_a);
% 
% % Enter a loop to brute force the solutions
% N = 1;
% for a=A
%     for b=B
%         p = [a b 0 0];
%         % Calculate the value
%         val = compensation_minFuncForBruteForce(p);
%         V(N) = val;
%         if mod(N,1000)==0
%             fprintf('N: %d\n',N)
%         end
%         N = N + 1;
%     end
% end
    
    
    










% fprintf('------> 1D\n')
% % Reserve arrays for storage
% min_a = 10000;
% max_a = 100000;
% n_a   = 1000;
% A = linspace(min_a,max_a,n_a);
% V = zeros(1,n_a);
% 
% % Enter a loop to brute force the solutions
% N = 1;
% for a=A
%     poly = [a 0 0];
%     % Calculate the value
%     val = compensation_minFuncForBruteForce_no_global(poly,w,I,p,Int_F);
%     V(N) = val;
%     if mod(N,1000)==0
%         fprintf('N: %d\n',N)
%     end
%     N = N + 1;
% end
% 
% figure(figNum)
%     figNum = figNum + 1;
%     plot(V)
    
    
    
    