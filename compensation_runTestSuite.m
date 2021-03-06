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

%% TEST 0: Starting conditions
fprintf('\n')
fprintf('TEST 0: Starting conditions\n')

w_low  = 1.82;
w_high = 1.84;
w_inc  = 0.0001;
w  = w_low:w_inc:w_high;
l  = 1e9*    2*pi*c ./ (w*1e15);
l  = fliplr(l);
p  = 80000*(w-w0).^2;
I  = exp(-(w-w0).^2/1e-5); %Caution: Probably the factor 1/lambda^2 is missing --> check

[Int_F,t_F,Ek_F] = compensation_calcFourierlimit(I,l);
[Int,t,E]        = compensation_calcFourierlimit(I,l,p);
factor = Int_F/Int;

figure(figNum)
    figNum = figNum + 1;
    [AX,H1,H2] = plotyy(w,I,w,p);
    set(AX,'xlim',[w_low w_high]);
    set(get(AX(1),'Ylabel'),'String','arbitrary units','fontweight','bold','fontsize',16)
    set(get(AX(2),'Ylabel'),'String','[rad]','fontweight','bold','fontsize',16)
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    title('Test 0: The pulse with intensity and spectrum we work with','fontweight','bold','fontsize',16);

%% TEST 1: See how our Fourier transformation performs.
fprintf('\n\nTEST 1: Fourier Trafo\n')
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
    xlabel('time [fs]','fontweight','bold','fontsize',16);
    ylabel('relative units','fontweight','bold','fontsize',16);
    legend('Fourier limit','Real pulse GDD 5e4','Real pulse GDD 8e4','Real pulse GDD 11e4')
    title('Test 1: Observation of the effects of different GDDs on the pulse','fontweight','bold','fontsize',16);
    
%% TEST 2: See if adding zeros increases our sample rate.
fprintf('\n\nTEST 2: Increase sample rate\n')

N = length(w);
[w_ext,I_ext,p_ext,l_ext] = compensation_extendPhaseByZeros(w,I,p,N);
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
    xlabel('time [fs]','fontweight','bold','fontsize',16);
    ylabel('relative units','fontweight','bold','fontsize',16);
    legend('Fourier limit 1xL','Real pulse 1xL','Fourier limit 2xL','Real pulse 2xL','Fourier limit 3xL','Real pulse 3xL')
    title('Test 2: Observation of the increase in peak sampling rate by adding zeros','fontweight','bold','fontsize',16);
    
%% TEST 3: See if our brute force ansatz works with our simple sample phase
fprintf('\n\nTEST 3: Approximation of the phase curve via brute force\n')

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
% CAREFUL: OUR POLYNOMIAL IS 1e4(w^2 - 2ww0 + w0^2) --> [8e4 -3.6576*8e4 3.3445*8e4] = [a b c]
a = 80000;
b = -2*w0*8e4;
c = w0^2*8e4;
for a=A
    polynomial = [a b c];
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
    xlabel('coefficient a [10^4]','fontweight','bold','fontsize',16);
    ylabel('Achieved minimum value','fontweight','bold','fontsize',16);
    title('Test 3: Calcultated values of our minimization function 1D case','fontweight','bold','fontsize',16);

fprintf('      > 2D\n')

% Enter a loop to brute force the solutions
p  = 1e6*(w-w0).^3+8e4*(w-w0).^2;
fprintf('p = 1e6*(w-w0).^3+8e4*(w-w0).^2;\n')
N = 1;

figure(figNum)
    figNum = figNum + 1;
    [AX,H1,H2] = plotyy(w,I,w,p);
    set(AX,'xlim',[w_low w_high]);
    set(get(AX(1),'Ylabel'),'String','arbitrary units','fontweight','bold','fontsize',16)
    set(get(AX(2),'Ylabel'),'String','[rad]','fontweight','bold','fontsize',16)
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    title('Test 3: New phase with p = 1e6*(w-w0).^3+8e4*(w-w0).^2','fontweight','bold','fontsize',16);

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
    h = surf(A*6,B*2,V); % The factors inserted here are our Taylor coefficients
    % Necessary, or the lines of the grid will overlay the colors
    set(h,'LineStyle','none')
    xlabel('TOD [fs^3]','fontweight','bold','fontsize',16);
    ylabel('GDD [fs^2]','fontweight','bold','fontsize',16);
    zlabel('Achieved minimum value (1-max(abs(E))','fontweight','bold','fontsize',16);
    title('Test 3: Calculated values of our minimization function 2D case','fontweight','bold','fontsize',16);

%% TEST 4: Apply the results on the measured phase
fprintf('\n\nTEST 4: Apply the results on the measured phase\n')

[t_Et,I_Et,p_Et,l_Sk,I_Sk,p_Sk]       = compensation_loadData('../Daten/Chirped Mirrors/Frogs/25.3W_RTT=6.404us_Ip=12.8A/','001_AIR_FROG_25.3W_RTT=6.404us_Ip=12.8A.bin');

[Int_F,t_F,Ek_F]                      = compensation_calcFourierlimit(I_Sk,l_Sk);

[l0,w0,w_spacing,w_Sk,I_Sk,p_Sk,l_Sk] = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,p_Sk,1030);

[B,A]                                 = compensation_makeFilter(w_Sk);

filtered_p_Sk                         = filtfilt(B,A,p_Sk);

[fit_w_Sk,fit_p_Sk,fit_I_Sk,fit_l_Sk] = compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,p_Sk);

[D1,D2,D3]                            = compensation_calcDevs(filtered_p_Sk,w_Sk,B,A);


figure(figNum)
    figNum = figNum + 1;
    [AX,H1,H2] = plotyy(w_Sk,I_Sk,w_Sk,p_Sk);
    set(AX,'xlim',[w_Sk(1) w_Sk(end)]);
    legend([H1;H2],'Spectrum','Phase')
    set(get(AX(1),'Ylabel'),'String','arbitrary units','fontweight','bold','fontsize',16)
    set(get(AX(2),'Ylabel'),'String','[rad]','fontweight','bold','fontsize',16)
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    title('Test 4: The pulse with intensity and spectrum as we measured it','fontweight','bold','fontsize',16);
    
figure(figNum)
    figNum = figNum + 1;
    plot(w_Sk,p_Sk)
    hold on
    plot(fit_w_Sk,fit_p_Sk,'LineWidth',2,'Color','red','LineStyle','-.')
    plot(w_Sk,I_Sk*40-20,'g')
    hold off
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    legend('Original','Fit range','Spectrum')
    title('Test 4: Visualization of the range we use to find our fit polynomial','fontweight','bold','fontsize',16);

polynomial = polyfit(fit_w_Sk,fit_p_Sk,4);
O4_fit     = polyval(polynomial,fit_w_Sk);
    
figure(figNum)
    figNum = figNum + 1;
    plot(fit_w_Sk,fit_p_Sk,'b')
    hold on
    plot(fit_w_Sk,fit_I_Sk*40-20,'g')
    plot(fit_w_Sk,O4_fit,'r')
    hold off
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    legend('Original','Spectrum','O(3) fit')
    title('Test 4: The resulting least square fit polynomial of degree 4','fontweight','bold','fontsize',16);
    
% Find the correct coefficients
C = polyval(polynomial,w0);
GD = polyval(polyder(polynomial),w0);
GDD = polyval(polyder(polyder(polynomial)),w0) / 2;
TOD = polyval(polyder(polyder(polyder(polynomial))),w0) / 6;
FOD = polyval(polyder(polyder(polyder(polyder(polynomial)))),w0) / 24;
fprintf(['\nAfter applying the following formulas:\n' ...
    'C = polyval(polynomial,w0)\n' ...
    'GD = polyval(polyder(polynomial),w0)\n' ... 
    'GDD = polyval(polyder(polyder(polynomial)),w0)\n' ...
    'TOD = polyval(polyder(polyder(polyder(polynomial))),w0)\n' ...
    'FOD = polyval(polyder(polyder(polyder(polyder(polynomial)))),w0)\n'])
fprintf('\nFrom our least square fit we get the following values:\n')
fprintf('FOD: %2.2e\n',FOD*24)
fprintf('TOD: %2.2e\n',TOD*6)
fprintf('GDD: %2.2e\n',GDD*2)
fprintf('GD: %2.2e\n',GD)
fprintf('C: %2.2e\n',C)

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
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    title('Test 4: View on the function bulk we run through compared to the original','fontweight','bold','fontsize',16);
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
    legend([h1 h2 h3],'Original phase minus the brute force solutions','Zero line','Original phase minus the least square fit')
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    title('Test 4: View on the differences we get when subtracting the estimate from the original','fontweight','bold','fontsize',16);
    hold off
    
N = 1;
for TOD2=A
    for GDD2=B
        polyn = [FOD TOD2 GDD2 GD C];
        % Calculate the value
        val = compensation_minFuncForBruteForce_no_global2(polyn,w0,w_Sk,I_Sk,p_Sk,Int_F);
        V(N) = val;
        if mod(N,1000)==0
            fprintf('N: %d\n',N)
        end
        N = N + 1;
    end
end

figure(figNum)
    figNum = figNum + 1;
    h = surf(A*6,B*2,V);
    % Necessary, or the lines of the grid will overlay the colors
    set(h,'LineStyle','none')
    xlabel('TOD [fs^3]','fontweight','bold','fontsize',16);
    ylabel('GDD [fs^2]','fontweight','bold','fontsize',16);
    zlabel('Achieved minimum value (1-max(abs(E))','fontweight','bold','fontsize',16);
    title('Test 4: Calculated 2D values at the vicinity of our least square fit','fontweight','bold','fontsize',16);
    
%% TEST 5: Comparing results with previous results from older functions
fprintf('\n\nTEST 5: Varify with results from older results\n')

C_main = -2.4188;
GD_main = 1947;
GDD_main = 41109;
TOD_main = -8.374e+05;
FOD_main = -1.8093e+07;

phase_main = polyval([FOD_main TOD_main GDD_main GD_main C_main],w_Sk-w0);

fprintf(['From our old function we found:\n' ...
    'C_main = %2.4e\n' ...
    'GD_main = %2.4e\n' ...
    'GDD_main = %2.4e\n' ...
    'TOD_main = %2.4e\n' ...
    'FOD_main = %2.4e\n'],C_main,GD_main,GDD_main*2,TOD_main*6,FOD_main*24)

figure(figNum)
    figNum = figNum + 1;
    plot(w_Sk,p_Sk)
    hold on
    plot(w_Sk,phase_main,'r')
    hold off
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    legend('Original','Old found optimum')
    title('Test 5: The resulting optimum fit polynomial of degree 4 from my old function','fontweight','bold','fontsize',16);
    
figure(figNum)
    figNum = figNum + 1;
    plot(w_Sk,p_Sk-phase_main)
    hold on
    plot(w_Sk,I_Sk*125,'g')
    hold off
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    legend('Difference','Spectrum')
    title('Test 5: The resulting difference when we subtract the old optimum from the original','fontweight','bold','fontsize',16);

% Append more zeros to increase the sample rate
N = 5*length(w_Sk);
[w_ext,I_ext,p_ext,l_ext] = compensation_extendPhaseByZeros(w_Sk,I_Sk,p_Sk-phase_main,N);
[Int_F,t_F,Ek_F] = compensation_calcFourierlimit(I_ext,l_ext);
[Int_ext,t_ext,E_ext] = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
factor = Int_F/Int_ext;
E_ext = E_ext * sqrt(factor);

figure(figNum)
    figNum = figNum + 1;
    plot(t_F,abs(Ek_F).^2,'g')
    hold on
    plot(t_ext,abs(E_ext).^2,'r')
    hold off
    xlim([-2000 2000])
    xlabel('time [fs]','fontweight','bold','fontsize',16);
    ylabel('relative units','fontweight','bold','fontsize',16);
    legend('Fourier limit','Compressed from old optimum')
    title('Test 5: Observation of the quality of compression on the pulse','fontweight','bold','fontsize',16);
    
val = compensation_minFuncForBruteForce_no_global2([FOD_main TOD_main GDD_main GD_main C_main],w0,w_Sk,I_Sk,p_Sk,Int_F);
fprintf('\nFound compression value: %2.2f%% of fourier limit\n',(1-val)*100)

%% TEST 6: Running the brute force attack with the knowledge from our old algorithm
fprintf('\n\nTEST 6: Finding a better brute force range with the knowledge from our old algorithm\n')

FOD = -1.8093e+07;
TOD = -8e+05;
GDD = 4e04;
GD = 1947;
C = -2.4188;

ord_TOD = floor(log10(abs(TOD)));
min_TOD = TOD-20*10^(floor(log10(abs(TOD)))-1);
max_TOD = TOD+20*10^(floor(log10(abs(TOD)))-1);
n_a   = 200;
ord_GDD = floor(log10(abs(GDD)));
min_GDD = GDD-5*10^(floor(log10(abs(GDD)))-1);
max_GDD = GDD+5*10^(floor(log10(abs(GDD)))-1);
n_b   = 50;
A = linspace(min_TOD,max_TOD,n_a);
B = linspace(min_GDD,max_GDD,n_b);
V = zeros(n_b,n_a);

fprintf('\nWe start with the following values:\n')
fprintf('FOD: %2.2e\n',FOD*24)
fprintf('TOD: [%2.2e  %2.2e]\n',min_TOD*6,max_TOD*6)
fprintf('GDD: [%2.2e  %2.2e]\n',min_GDD*2,max_GDD*2)
fprintf('GD: %2.2e\n',GD)
fprintf('C: %2.2e\n',C)

fprintf('Using the following calculations:\n')
fprintf(['ord_TOD = floor(log10(abs(TOD)))\n' ...
'min_TOD = TOD-20*10^(floor(log10(abs(TOD)))-1)\n' ...
'max_TOD = TOD+20*10^(floor(log10(abs(TOD)))-1)\n' ...
'n_a   = 200\n' ...
'ord_GDD = floor(log10(abs(GDD)))\n' ...
'min_GDD = GDD-5*10^(floor(log10(abs(GDD)))-1)\n' ...
'max_GDD = GDD+5*10^(floor(log10(abs(GDD)))-1)\n' ...
'n_b   = 50\n' ...
'A = linspace(min_TOD,max_TOD,n_a)\n' ...
'B = linspace(min_GDD,max_GDD,n_b)\n' ...
'V = zeros(n_b,n_a)\n'])

figure(figNum)
    figNum = figNum + 1;
    hold on
for TOD2=A
    for GDD2=B
        polyn = [FOD TOD2 GDD2 GD C];
        plot(w_Sk,polyval(polyn,w_Sk-w0))
    end
end
    h1 = plot(w_Sk,polyval(polyn,w_Sk-w0));
    h2 = plot(w_Sk,p_Sk,'r');
    legend([h1 h2],'Brute force','Original')
%     ylim([min(fit_p_Sk) max(fit_p_Sk)])
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    title('Test 6: View on the function bulk we run through compared to the original','fontweight','bold','fontsize',16);
    hold off
    
figure(figNum)
    figNum = figNum + 1;
    hold on
for TOD2=A
    for GDD2=B
        polyn = [FOD TOD2 GDD2 GD C];
        plot(w_Sk,p_Sk-polyval(polyn,w_Sk-w0))
    end
end
    h1 = plot(w_Sk,p_Sk-polyval(polyn,w_Sk-w0));
    h2 = plot(w_Sk,zeros(size(w_Sk)),'r');
    h3 = plot(w_Sk,p_Sk-polyval([FOD TOD GDD GD C],w_Sk-w0),'g');
    legend([h1 h2 h3],'Original phase minus the brute force solutions','Zero line','Original phase minus the least square fit')
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    title('Test 6: View on the differences we get when subtracting the estimate from the original','fontweight','bold','fontsize',16);
    hold off
    
N = 1;
for TOD2=A
    for GDD2=B
        polyn = [FOD TOD2 GDD2 GD C];
        % Calculate the value
        val = compensation_minFuncForBruteForce_no_global2(polyn,w0,w_Sk,I_Sk,p_Sk,Int_F);
        V(N) = val;
        if mod(N,1000)==0
            fprintf('N: %d\n',N)
        end
        N = N + 1;
    end
end

figure(figNum)
    figNum = figNum + 1;
    h = surf(A*6,B*2,V);
    % Necessary, or the lines of the grid will overlay the colors
    set(h,'LineStyle','none')
    xlabel('TOD [fs^3]','fontweight','bold','fontsize',16);
    ylabel('GDD [fs^2]','fontweight','bold','fontsize',16);
    zlabel('Achieved minimum value (1-max(abs(E))','fontweight','bold','fontsize',16);
    title('Test 6: Calculated values of our minimization function 2D case','fontweight','bold','fontsize',16);
    
%% TEST 7: Search for an optimal solution
fprintf('\n\nTEST 7: Search for an optimal solution using optimization algorithms\n')

fprintf('\nWe start with the following values:\n')
fprintf('FOD: %2.2e\n',FOD*24)
fprintf('TOD: %2.2e\n',TOD*6)
fprintf('GDD: %2.2e\n',GDD*2)
fprintf('GD: %2.2e\n',GD)
fprintf('C: %2.2e\n',C)

options = optimset('MaxIter', 1000,'MaxFunEvals',1e5);
min_func = @(x)compensation_minFuncForBruteForce_no_global2([FOD x(1) x(2) GD C],w0,w_Sk,I_Sk,p_Sk,Int_F);

[solution,val] = fminsearch(min_func,[TOD GDD],options);

fprintf('With the new optimization we find a value of %2.2f%% of the fourier peak\n',(1-val)*100)
fprintf('The found values for TOD and GDD are:\n')
fprintf('TOD: %2.2e\n',solution(1)*6)
fprintf('GDD: %2.2e\n',solution(2)*2)

phase_new = polyval([FOD solution(1) solution(2) GD C],w_Sk-w0);
N = 5*length(w_Sk);
[w_ext,I_ext,p_ext,l_ext] = compensation_extendPhaseByZeros(w_Sk,I_Sk,p_Sk-phase_new,N);
[Int_F,t_F,Ek_F] = compensation_calcFourierlimit(I_ext,l_ext);
[Int_new,t_new,E_new] = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
factor = Int_F/Int_new;

figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(w_Sk,p_Sk,'r');
    plot(w_Sk,phase_new,'b');
    legend('Original','New found optimum')
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    title('Test 7: Comparison of the new found optimum to the original phase','fontweight','bold','fontsize',16);
    hold off
    
figure(figNum)
    figNum = figNum + 1;
    plot(w_Sk,p_Sk-phase_new);
    legend('Original phase minus new found optimum')
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    title('Test 7: Difference between original phase and new optimum','fontweight','bold','fontsize',16);

figure(figNum)
    figNum = figNum + 1;
    plot(t_F,abs(Ek_F).^2,'g')
    hold on
    plot(t_new,abs(E_new).^2.*factor,'r')
    plot(t_ext,abs(E_ext).^2,'b')
    hold off
    xlim([-2000 2000])
    xlabel('time [fs]','fontweight','bold','fontsize',16);
    ylabel('relative units','fontweight','bold','fontsize',16);
    legend('Fourier limit','Compressed from new optimum','Compressed from old optimum')
    title('Test 7: Observation of the compression of our new found optimum','fontweight','bold','fontsize',16);
    
%% TEST 8: Considering higher order optimization
fprintf('\n\nTEST 8: Including higher orders in the optimization process\n')

fprintf('\nWe start again with the following values:\n')
fprintf('FOD: %2.2e\n',FOD*24)
fprintf('TOD: %2.2e\n',TOD*6)
fprintf('GDD: %2.2e\n',GDD*2)
fprintf('GD: %2.2e\n',GD)
fprintf('C: %2.2e\n',C)

options = optimset('MaxIter', 1000,'MaxFunEvals',1e5);
min_func = @(x)compensation_minFuncForBruteForce_no_global2([x(1) x(2) x(3) GD C],w0,w_Sk,I_Sk,p_Sk,Int_F);

[solution,val] = fminsearch(min_func,[FOD TOD GDD],options);

fprintf('Including the 4th order in the optimization we find a value of %2.2f%% of the fourier peak\n',(1-val)*100)
fprintf('The found values for FOD, TOD and GDD are:\n')
fprintf('FOD: %2.2e\n',solution(1)*24)
fprintf('TOD: %2.2e\n',solution(2)*6)
fprintf('GDD: %2.2e\n',solution(3)*2)

phase_new = polyval([solution(1) solution(2) solution(3) GD C],w_Sk-w0);
N = 5*length(w_Sk);
[w_ext,I_ext,p_ext,l_ext] = compensation_extendPhaseByZeros(w_Sk,I_Sk,p_Sk-phase_new,N);
[Int_F,t_F,Ek_F] = compensation_calcFourierlimit(I_ext,l_ext);
[Int_new2,t_new2,E_new2] = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
factor2 = Int_F/Int_new2;
% Consider only the narrow spectral range of interest
narrow_phase_new = polyval([solution(1) solution(2) solution(3) GD C],fit_w_Sk-w0);
N = 5*length(fit_w_Sk);
[narrow_w_ext,narrow_I_ext,narrow_p_ext,narrow_l_ext] = compensation_extendPhaseByZeros(fit_w_Sk,fit_I_Sk,fit_p_Sk-narrow_phase_new,N);
[narrow_Int_F,narrow_t_F,narrow_Ek_F] = compensation_calcFourierlimit(narrow_I_ext,narrow_l_ext);
[narrow_Int_new2,narrow_t_new2,narrow_E_new2] = compensation_calcFourierlimit(narrow_I_ext,narrow_l_ext,narrow_p_ext);
narrow_factor2 = narrow_Int_F/narrow_Int_new2;
fprintf('Restricting the calculations to the narrow spectral range leads %2.2f%% of the fourier peak\n',narrow_factor2*100)

figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(w_Sk,p_Sk,'r');
    plot(w_Sk,phase_new,'b');
    legend('Original','New found optimum')
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    title('Test 8: Comparison of the found optimum including 4th order to the original phase','fontweight','bold','fontsize',16);
    hold off
    
figure(figNum)
    figNum = figNum + 1;
    plot(w_Sk,p_Sk-phase_new);
    legend('Original phase minus new found optimum')
    xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
    ylabel('[rad]','fontweight','bold','fontsize',16);
    title('Test 8: Difference between original phase and compression including order 4','fontweight','bold','fontsize',16);

figure(figNum)
    figNum = figNum + 1;
    plot(t_F,abs(Ek_F).^2,'g')
    hold on
    plot(t_new,abs(E_new).^2.*factor,'r')
    plot(t_new2,abs(E_new2).^2.*factor2,'c')
    plot(t_ext,abs(E_ext).^2,'b')
    hold off
    xlim([-2000 2000])
    xlabel('time [fs]','fontweight','bold','fontsize',16);
    ylabel('relative units','fontweight','bold','fontsize',16);
    legend('Fourier limit','Compression including up to order 3','Compression including up to order 4','Old found optimum')
    title('Test 8: Observation of the compression with higher orders','fontweight','bold','fontsize',16);
    
%% TEST 9: Showing the effect of higher orders
fprintf('\n\nTEST 9: Investigating the effect of higher orders to the final compression\n')

orders = 2:10;
vals = zeros(1,length(orders));
options = optimset('MaxIter', 1000,'MaxFunEvals',1e5);
min_func = @(x)compensation_minFuncForBruteForce_no_global2(x,w0,w_Sk,I_Sk,p_Sk,Int_F);

for order=orders
    polynomial = polyfit(fit_w_Sk,fit_p_Sk,order);
    coeffs = zeros(1,order+1);

    coeffs(end) = polyval(polynomial,w0);
    deriv = polyder(polynomial);
    coeffs(end-1) = polyval(deriv,w0);
    for i=2:order
        deriv = polyder(deriv);
        coeffs(end-i) = polyval(deriv,w0) / factorial(i);
    end
    fprintf('Order %d\n',order)
    fprintf('Starting values:')
    disp(coeffs)
    [solution,val] = fminsearch(min_func,coeffs,options);
    fprintf('Final values:')
    disp(solution)
    vals(order-1) = 1-val;
end

figure(figNum)
    figNum = figNum + 1;
    plot(orders,vals)
    xlabel('Highest included order','fontweight','bold','fontsize',16);
    ylabel('Achieved percentage of fourier limit','fontweight','bold','fontsize',16);
    title('Test 9: Observation of the effectiveness to include orders higher than 4','fontweight','bold','fontsize',16);