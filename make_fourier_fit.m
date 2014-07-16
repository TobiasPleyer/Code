function [solution,val]=make_fourier_fit(omega,intensity,phase,fourierlimit,poly)
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
    c = 299792458;
    lambda = 2*pi*c ./ (omega*1e15);
%     p = polyfit(omega,phase,order)
    min_func = @(x)observe_compensation_min_Func(x,omega,intensity,phase,lambda',fourierlimit);
    [solution,val] = fminsearch(min_func,poly);
    % Now we have to adjust the linear and constant term so that we can see
    % a comparison with the other phases
    min_func = @(x)sum((polyval(poly,omega)-polyval(solution,omega)-polyval(x,omega)).^2);
    x0 = [0,0];
    [s,v] = fminsearch(min_func,x0);
    fprintf('Optimum value from make_fourier_fit: %2.2f\n',v)
    disp(s)
    solution(end-1:end) = solution(end-1:end) + s;
end