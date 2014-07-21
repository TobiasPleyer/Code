function [solution,val]=make_fourier_fit(omega,intensity,phase,fourierlimit,poly,figNum)
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
    
    fval_history = [];
    iter_history = [];
    
    c = 299792458;
    lambda = 2*pi*c ./ (omega*1e15);

    options = optimset('MaxIter', 1000,'MaxFunEvals',1e5,'OutputFcn',@outfunc);
    min_func = @(x)observe_compensation_min_Func(x,omega,intensity,phase,lambda',fourierlimit);
    
    [solution,val] = fminsearch(min_func,poly(1:end-2),options);
    fprintf('For order %d we found: %2.2f with start values (%2.3e,%2.3e)\n',length(poly)-1,1-val,poly(1),poly(2))
    fprintf('Final constants:\n')
    disp(solution)
    
    % Don't forget to append back the last two constants at the end
    solution = [solution poly(end-1:end)];
    
    % Now we have to adjust the linear and constant term so that we can see
    % a comparison with the other phases
    min_func = @(x)sum((polyval(poly,omega)-polyval(solution,omega)-polyval(x,omega)).^2);
    x0 = [0,0];
    [s,v] = fminsearch(min_func,x0);
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
 
    figure(figNum)
    plot(iter_history,fval_history,'x')
    xlabel('Iteration #')
    ylabel('Returned optimum')
    title(sprintf('Monitoring of the optimization trend curve for order %d.',length(poly)-1))
    
    % Call test_initial_conditions.m
    test_initial_conditions(omega,lambda,intensity,phase,fourierlimit,poly,figNum+1,solution)
end