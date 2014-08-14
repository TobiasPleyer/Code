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
GD  = 299792458;
l0 = 1030;
w0 = 1e-15*  2*pi*GD / (l0*1e-9);

% TEST 0: Starting conditions
fprintf('\n')
fprintf('TEST 0: Starting conditions\n')

w_low  = 1.82;
w_high = 1.84;
w_inc  = 0.0001;
w  = w_low:w_inc:w_high;
l  = 1e9*    2*pi*GD ./ (w*1e15);
l  = fliplr(l);
p  = 80000*(w-w0).^2;
I  = exp(-(w-w0).^2/1e-5);

[Int_F,t_F,Ek_F] = compensation_calcFourierlimit(I,l);
[Int,t,E]        = compensation_calcFourierlimit(I,l,p);
factor = Int_F/Int;

figure(figNum)
    figNum = figNum + 1;
    [AX,H1,H2] = plotyy(w,I,w,p);
    set(AX,'xlim',[w_low w_high]);
    set(get(AX(1),'Ylabel'),'String','arbitrary units')
    set(get(AX(2),'Ylabel'),'String','[rad]')
    xlabel('omega [1/fs]')
    title('Test 0: The pulse with intensity and spectrum we work with')

% TEST 1: See how our Fourier transformation performs.
fprintf('TEST 1: Fourier Trafo\n')
fprintf('p1 = 50000*(w-w0).^2\np2 = 80000*(w-w0).^2\np3 = 110000*(w-w0).^2;\n')

p1  = 50000*(w-w0).^2;
p2  = 80000*(w-w0).^2;
p3  = 110000*(w-w0).^2;
[Int1,t1,E1] = compensation_calcFourierlimit(I,l,p1);
[Int2,t2,E2] = compensation_calcFourierlimit(I,l,p2);
[Int3,t3,E3] = compensation_calcFourierlimit(I,l,p3);
factor1 = Int_F/Int1;
factor2 = Int_F/Int2;
factor3 = Int_F/Int3;

figure(figNum)
    figNum = figNum + 1;
    plot(t_F,abs(Ek_F).^2,'g')
    hold on
    plot(t1,abs(E1).^2.*factor1,'r')
    plot(t2,abs(E2).^2.*factor2,'b')
    plot(t3,abs(E3).^2.*factor3,'m')
    hold off
    xlim([-2000 2000])
    xlabel('time [fs]')
    ylabel('relative units')
    legend('Fourier limit','Real pulse GDD 5e4','Real pulse GDD 8e4','Real pulse GDD 11e4')
    title('Test 1: Observation of the effects of different GDDs on the pulse')
    
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
    title('Test 2: Observation of the increase in peak sampling rate by adding zeros')
    
% TEST 3: See if our brute force ansatz works with our simple sample phase
fprintf('TEST 3: Approximation of the phase curve via brute force\n')

fprintf('      > 1D\n')
p = 80000*(w-w0).^2;
% Reserve arrays for storage
min_TOD = 10000;
max_TOD = 100000;
n_a   = 10;
A = linspace(min_TOD,max_TOD,n_a);
V = zeros(1,n_a);

% Enter a loop to brute force the solutions
N = 1;
% CAREFUL: OUR POLYNOMIAL IS 1e4(w^2 - 2ww0 + w0^2) --> [8e4 -3.6576*8e4 3.3445*8e4]
TOD = 80000;
GDD = -2*w0*8e4;
GD = w0^2*8e4;
for TOD=A
    polynomial = [TOD GDD GD];
    % Calculate the value
    val = compensation_minFuncForBruteForce_no_global(polynomial,w,I,p,Int_F);
    V(N) = val;
    if mod(N,1000)==0
        fprintf('N: %d\n',N)
    end
    N = N + 1;
end

figure(figNum)
    figNum = figNum + 1;
    plot(V)
    title('Test 3: Calcultated values of our minimization function 1D case')

fprintf('      > 2D\n')

% Enter a loop to brute force the solutions
p  = 1e6*(w-w0).^3+8e4*(w-w0).^2;
N = 1;

figure(figNum)
    figNum = figNum + 1;
    [AX,H1,H2] = plotyy(w,I,w,p);
    set(AX,'xlim',[w_low w_high]);
    set(get(AX(1),'Ylabel'),'String','arbitrary units')
    set(get(AX(2),'Ylabel'),'String','[rad]')
    xlabel('omega [1/fs]')
    title('Test 3: New phase with p = 5000*(w-w0).^3+80000*(w-w0).^2')

% Reserve arrays for storage
min_TOD = 5e5;
max_TOD = 5e6;
n_a   = 50;
min_GDD = 1e4;
max_GDD = 1e5;
n_b   = 50;
A = linspace(min_TOD,max_TOD,n_a);
B = linspace(min_GDD,max_GDD,n_b);
V = zeros(n_b,n_a);

for TOD=A
    for GDD=B
        polynomial = compensation_calcPoly3([TOD GDD],w0);
        % Calculate the value
        val = compensation_minFuncForBruteForce_no_global(polynomial,w,I,p,Int_F);
        V(N) = val;
%         if mod(N,1000)==0
%             fprintf('N: %d\n',N)
%         end
        N = N + 1;
    end
end

figure(figNum)
    figNum = figNum + 1;
    h = surf(A,B,V);
    % Necessary, or the lines of the grid will overlay the colors
    set(h,'LineStyle','none')
    title('Test 3: Calcultated values of our minimization function 2D case')

% TEST 4: Apply the results on the measured phase
fprintf('TEST 4: Apply the results on the measured phase\n')

[t_Et,I_Et,p_Et,l_Sk,I_Sk,p_Sk]       = compensation_loadData();

[Int_F,t_F,Ek_F]                      = compensation_calcFourierlimit(I_Sk,l_Sk);

[l0,w0,w_spacing,w_Sk,I_Sk,p_Sk,l_Sk] = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,p_Sk);

[B,A]                                 = compensation_makeFilter(w_Sk);

filtered_p_Sk                         = filtfilt(B,A,p_Sk);

[fit_w_Sk,fit_p_Sk,fit_I_Sk,fit_l_Sk] = compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,p_Sk);

[D1,D2,D3]                            = compensation_calcDevs(filtered_p_Sk,w_Sk,B,A);


figure(figNum)
    figNum = figNum + 1;
    [AX,H1,H2] = plotyy(w_Sk,I_Sk,w_Sk,p_Sk);
    set(AX,'xlim',[w_Sk(1) w_Sk(end)]);
    legend([H1;H2],'Spectrum','Phase')
    set(get(AX(1),'Ylabel'),'String','arbitrary units')
    set(get(AX(2),'Ylabel'),'String','[rad]')
    xlabel('omega [1/fs]')
    title('Test 4: The pulse with intensity and spectrum as we measured it')
    
figure(figNum)
    figNum = figNum + 1;
    plot(w_Sk,p_Sk)
    hold on
    plot(fit_w_Sk,fit_p_Sk,'LineWidth',2,'Color','red','LineStyle','-.')
    plot(w_Sk,I_Sk*40-20,'g')
    hold off
    xlabel('omega [1/fs]')
    ylabel('[rad]')
    legend('Original','Fit range','Spectrum')
    title('Test 4: Visualization of the range we use to find our fit polynomial')

polynomial = polyfit(fit_w_Sk,fit_p_Sk,4);
O4_fit     = polyval(polynomial,fit_w_Sk);
    
figure(figNum)
    figNum = figNum + 1;
    plot(fit_w_Sk,fit_p_Sk,'b')
    hold on
    plot(fit_w_Sk,fit_I_Sk*40-20,'g')
    plot(fit_w_Sk,O4_fit,'r')
    hold off
    xlabel('omega [1/fs]')
    ylabel('[rad]')
    legend('Original','Spectrum','O(3) fit')
    title('Test 4: The resulting least square fit polynomial of degree 4')
    
% Find the correct coefficients
C = polyval(polynomial,w0);
GD = polyval(polyder(polynomial),w0);
GDD = polyval(polyder(polyder(polynomial)),w0) / 2;
TOD = polyval(polyder(polyder(polyder(polynomial))),w0) / 6;
FOD = polyval(polyder(polyder(polyder(polyder(polynomial)))),w0) / 24;

% One option
% ord_a = floor(log10(abs(a)));
% n_a   = 50;
% ord_b = floor(log10(abs(b)));
% n_b   = 50;
% A = sign(a).*logspace(ord_a-1,ord_a+1,n_a);
% B = sign(b).*logspace(ord_b-1,ord_b+1,n_b);
% V = zeros(n_b,n_a);

ord_TOD = floor(log10(abs(TOD)));
min_TOD = TOD-1*10^(floor(log10(abs(TOD)))-1);
max_TOD = TOD+1*10^(floor(log10(abs(TOD)))-1);
n_a   = 50;
ord_GDD = floor(log10(abs(GDD)));
min_GDD = GDD-1*10^(floor(log10(abs(GDD)))-1);
max_GDD = GDD+1*10^(floor(log10(abs(GDD)))-1);
n_b   = 50;
A = linspace(min_TOD,max_TOD,n_a);
B = linspace(min_GDD,max_GDD,n_b);
V = zeros(n_b,n_a);

figure(figNum)
    figNum = figNum + 1;
    hold on
for TOD2=A
    for GDD2=B
        polyn = [FOD TOD2 GDD2 GD C];
        plot(fit_w_Sk,polyval(polyn,fit_w_Sk-w0))
    end
end
    h1 = plot(fit_w_Sk,polyval(polyn,fit_w_Sk-w0));
    h2 = plot(fit_w_Sk,fit_p_Sk,'r');
    legend([h1 h2],'Brute force','Original')
    ylim([min(fit_p_Sk) max(fit_p_Sk)])
    xlabel('omega [1/fs]')
    ylabel('[rad]')
    title('View on the function bulk we run through compared to the original')
    hold off
    
figure(figNum)
    figNum = figNum + 1;
    hold on
for TOD2=A
    for GDD2=B
        polyn = [FOD TOD2 GDD2 GD C];
        plot(fit_w_Sk,fit_p_Sk-polyval(polyn,fit_w_Sk-w0))
    end
end
    h1 = plot(fit_w_Sk,fit_p_Sk-polyval(polyn,fit_w_Sk-w0));
    h2 = plot(fit_w_Sk,zeros(size(fit_w_Sk)),'r');
    h3 = plot(fit_w_Sk,fit_p_Sk-polyval([FOD TOD GDD GD C],fit_w_Sk-w0),'g');
    legend([h1 h2 h3],'Brute force','Zero line','Least square')
    xlabel('omega [1/fs]')
    ylabel('[rad]')
    title('View on the differences we get when subtracting the estimate from the original')
    hold off
    
N = 1;
for TOD2=A
    for GDD2=B
        polyn = [FOD TOD2 GDD2 GD C];
        % Calculate the value
        val = compensation_minFuncForBruteForce_no_global2(polyn,w0,w,I,p,Int_F);
        V(N) = val;
        if mod(N,1000)==0
            fprintf('N: %d\n',N)
        end
        N = N + 1;
    end
end

figure(figNum)
    figNum = figNum + 1;
    h = surf(A,B,V);
    % Necessary, or the lines of the grid will overlay the colors
    set(h,'LineStyle','none')
    title('Test 4: Calcultated values of our minimization function 2D case')
    
% TEST 5: Comparing results with previous resluts from older functions
fprintf('TEST 5: Varify with results from older results\n')

C_main = -2.4188;
GD_main = 1947;
GDD_main = 41109;
TOD_main = -8.374e+05;
FOD_main = -1.8093e+07;

phase_main = polyval([FOD_main TOD_main GDD_main GD_main C_main],w_Sk-w0);

fprintf('From our old function we found:\nC_main = -2.4188\nGD_main = 1947\nGDD_main = 41109\nTOD_main = -8.374e+05\nFOD_main = -1.8093e+07\n')

figure(figNum)
    figNum = figNum + 1;
    plot(w_Sk,p_Sk)
    hold on
    plot(w_Sk,phase_main,'r')
    hold off
    xlabel('omega [1/fs]')
    ylabel('[rad]')
    legend('Original','Old found optimum')
    title('Test 5: The resulting optimum fit polynomial of degree 4 from my old function')
    
figure(figNum)
    figNum = figNum + 1;
    plot(w_Sk,p_Sk-phase_main)
    hold on
    plot(w_Sk,I_Sk*125,'g')
    hold off
    xlabel('omega [1/fs]')
    ylabel('[rad]')
    legend('Difference','Spectrum')
    title('Test 5: The resulting difference when we subtract the old optimum from the original')
    
[Int_main,t_main,E_main] = compensation_calcFourierlimit(I_Sk,fliplr(l_Sk'),p_Sk-phase_main);
factor = Int_F/Int_main;

figure(figNum)
    figNum = figNum + 1;
    plot(t_F,abs(Ek_F).^2,'g')
    hold on
    plot(t_main,abs(E_main).^2.*factor,'r')
    hold off
    xlim([-2000 2000])
    xlabel('time [fs]')
    ylabel('relative units')
    legend('Fourier limit','Compressed from old optimum')
    title('Test 5: Observation of the quality of compression on the pulse')
    
val = compensation_minFuncForBruteForce_no_global2([FOD_main TOD_main GDD_main GD_main C_main],w0,w_Sk,I_Sk,p_Sk,Int_F);
fprintf('\nFound compression value: %2.2f%% of fourier limit\n',(1-val)*100)
   
    





    
    
    
    