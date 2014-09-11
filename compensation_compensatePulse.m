function [w_Sk,I_Sk,p_Sk,phase_new]=compensation_compensatePulse(folder,fileBase,pm,varargin)

    if nargin < 3
        error('Usage: [w_Sk,I_Sk,p_Sk,phase_new]=compensation_compensatePulse(folder,fileBase,pm,varargin)')
    end
    
    if nargin == 3
        [~,~,~,l_Sk,I_Sk,p_Sk]       = compensation_loadData(folder,fileBase);
        [Int_F,~,~]                  = compensation_calcFourierlimit(I_Sk,l_Sk);
        [~,w0,~,w_Sk,I_Sk,p_Sk,l_Sk] = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,pm*p_Sk);
        [fit_w_Sk,fit_p_Sk,~,~]      = compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,p_Sk);
        polynomial                   = polyfit(fit_w_Sk,fit_p_Sk,4);
        % Find the correct coefficients
        C   = polyval(polynomial,w0);
        GD  = polyval(polyder(polynomial),w0);
        GDD = polyval(polyder(polyder(polynomial)),w0);
        TOD = polyval(polyder(polyder(polyder(polynomial))),w0);
    %     FOD = polyval(polyder(polyder(polyder(polyder(polynomial)))),w0);

        fprintf('\nWe start with the following values found from our least square approximation:\n')
        fprintf('TOD: %2.2e\n',TOD)
        fprintf('GDD: %2.2e\n',GDD)
        fprintf('GD: %2.2e\n',GD)
        fprintf('C: %2.2e\n',C)

        options = optimset('MaxIter', 1000,'MaxFunEvals',1e5);
        min_func = @(x)compensation_minFuncForBruteForce_no_global2([x(1) x(2) GD C],w0,w_Sk,I_Sk,p_Sk,Int_F);

        [solution,val] = fminsearch(min_func,[TOD/6 GDD/2],options);

        fprintf('\nWith our optimization routine we find a value of %2.2f%% of the fourier peak\n',(1-val)*100)
        fprintf('The found values for TOD and GDD are:\n')
        fprintf('TOD: %2.2e\n',solution(1)*6)
        fprintf('GDD: %2.2e\n',solution(2)*2)

        phase_new             = polyval([solution(1) solution(2) GD C],w_Sk-w0);
        N                     = 5*length(w_Sk);
        [~,I_ext,p_ext,l_ext] = compensation_extendPhaseByZeros(w_Sk,I_Sk,p_Sk-phase_new,N);
        [Int_F,t_F,Ek_F]      = compensation_calcFourierlimit(I_ext,l_ext);
        [Int_new,t_new,E_new] = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
        [~,I_ext,p_ext,l_ext] = compensation_extendPhaseByZeros(w_Sk,I_Sk,p_Sk,N);
        [Int_orig,t_orig,E_orig] = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
        factor                = Int_F/Int_new;
        factor2               = Int_F/Int_orig;

        figure()
            hold on
            plot(w_Sk,I_Sk*20,'g')
            plot(w_Sk,p_Sk,'r');
            plot(w_Sk,phase_new,'b');
            legend('Spectrum (intensity)','Original phase','Approximated phase')
            ylabel('[rad]','fontweight','bold','fontsize',16)
            xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
            ylim([-5 20])
            title('Comparison of the found optimum to the original phase','fontweight','bold','fontsize',16);
            hold off

        figure()
            plot(t_F,abs(Ek_F).^2,'g')
            hold on
            plot(t_new,abs(E_new).^2.*factor,'r')
            plot(t_orig,abs(E_orig).^2.*factor2,'b')
            hold off
            xlim([-2000 2000])
            xlabel('time [fs]','fontweight','bold','fontsize',16);
            ylabel('relative units','fontweight','bold','fontsize',16);
            legend('Fourier limit','Compressed from new optimum','Before compensation')
            title('Observation of the compression of our found optimum','fontweight','bold','fontsize',16);
    end
    
    if nargin == 4
        [~,~,~,l_Sk,I_Sk,p_Sk]      = compensation_loadData(folder,fileBase);
        [~,~,~]                     = compensation_calcFourierlimit(I_Sk,l_Sk);
        [~,~,~,w_Sk,I_Sk,p_Sk,l_Sk] = compensation_makeOmegaEqualSpaced(l_Sk,I_Sk,pm*p_Sk);
        [~,~,~,~]                   = compensation_findRangeOfInterest(I_Sk,w_Sk,l_Sk,p_Sk);
        N                           = 5*length(w_Sk);
        phase_new                   = varargin{1};
        [~,I_ext,p_ext,l_ext]       = compensation_extendPhaseByZeros(w_Sk,I_Sk,p_Sk-phase_new,N);
        [Int_F,t_F,Ek_F]            = compensation_calcFourierlimit(I_ext,l_ext);
        [Int_new,t_new,E_new]       = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
        [~,I_ext,p_ext,l_ext]       = compensation_extendPhaseByZeros(w_Sk,I_Sk,p_Sk,N);
        [Int_orig,t_orig,E_orig]    = compensation_calcFourierlimit(I_ext,l_ext,p_ext);
        factor                      = Int_F/Int_new;
        factor2                     = Int_F/Int_orig;

        figure()
            hold on
            plot(w_Sk,I_Sk*20,'g')
            plot(w_Sk,p_Sk,'r');
            plot(w_Sk,phase_new,'b');
            legend('Spectrum (intensity)','Original phase','Approximated phase')
            ylabel('[rad]','fontweight','bold','fontsize',16)
            xlabel('omega [1/fs]','fontweight','bold','fontsize',16);
            ylim([-5 20])
            title('Comparison of the found optimum to the original phase','fontweight','bold','fontsize',16);
            hold off

        figure()
            plot(t_F,abs(Ek_F).^2,'g')
            hold on
            plot(t_new,abs(E_new).^2.*factor,'r')
            plot(t_orig,abs(E_orig).^2.*factor2,'b')
            hold off
            xlim([-2000 2000])
            xlabel('time [fs]','fontweight','bold','fontsize',16);
            ylabel('relative units','fontweight','bold','fontsize',16);
            legend('Fourier limit','Compressed from new optimum','Before compensation')
            title('Observation of the compression of our found optimum','fontweight','bold','fontsize',16);
    end
end