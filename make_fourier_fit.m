function [solution,val]=make_fourier_fit()
    %% ########################################################################
    %%      Function make_fourier_fit
    %
    %       This function is a helper function to use with fminsearch.
    %       It takes the angular frequency, intensity and phase of a pulse 
    %       as well as its Fourier limit and return the difference 1-peak,
    %       were peak is normed with respect to the Fourier limit.
    %  
    %  Version: 1.0
    %
    %  INPUTS: [omega,intensity,phase,fourierlimit,order]
    %       omega:      Wavelength in 1/fs
    %       intensity:  Intensity in arbitrary units
    %       phase:      Phase in radians
    %       fourierlimit: The value of the integral of the Fourier limit
    %       order:      The order of the polynomial fitted as an integer
    %
    %  OUTPUTS: [solution,val]
    %       solution:   The coefficients of the found polynomial as Matlab
    %                   defines it (can be used with polyval etc...)
    %       val:        The value of the achieved optimization. 1-val gives
    %                   the achieved percentage value of the Fourier limit.
    %
    %  HISTORY:
    %      v1.0: First runable state. (15.07.2014)
    %
    %% ########################################################################
    
    %% Too make function definitions easier we define a bunch of global variables

    global w_Sk
    global p figNum
    
    %%


    %% Code
    
    fval_history = [];
    iter_history = [];

    options = optimset('MaxIter', 1000,'MaxFunEvals',1e5,'OutputFcn',@outfunc);
    min_func = @(x)observe_compensation_min_Func(x);
    
    [solution,val] = fminsearch(min_func,p(1:end-2),options);
    fprintf('For order %d we found: %2.2f with start values (%2.3e,%2.3e)\n',length(p)-1,1-val,p(1),p(2))
    fprintf('Final constants:\n')
    disp(solution)
    
    % Don't forget to append back the last two constants at the end
    solution = [solution p(end-1:end)];
    
    % Now we have to adjust the linear and constant term so that we can see
    % a comparison with the other phases
    min_func = @(x)sum((polyval(p,w_Sk)-polyval(solution,w_Sk)-polyval(x,w_Sk)).^2);
    x0 = [0,0];
    options = optimset('MaxIter', 1000,'MaxFunEvals',1e5,'OutputFcn',@outfunc2);
    [s,~] = fminsearch(min_func,x0,options);
    solution(end-1:end) = solution(end-1:end) + s;
    
    function stop=outfunc(x,optimvalues,state)
        stop = false;
        fval = optimvalues.fval;
        iteration = optimvalues.iteration;
        if isequal(state,'iter')
          fval_history = [fval_history; fval];
          iter_history = [iter_history; iteration];
        end
    end

    function stop=outfunc2(x,optimvalues,state)
        stop = false;
        fprintf('x: [%2.10f, %2.10f]\n',x(1),x(2))
    end
 
%     figure(figNum)
%     figNum = figNum + 1;
%     plot(iter_history,fval_history,'x')
%     xlabel('Iteration #')
%     ylabel('Returned optimum')
%     title(sprintf('Monitoring of the optimization trend curve for order %d.',length(p)-1))
    
    % Call test_initial_conditions.m
%     test_initial_conditions(solution)
end