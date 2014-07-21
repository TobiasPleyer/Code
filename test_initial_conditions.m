function test_initial_conditions(omega,lambda,intensity,phase,fourierlimit,poly,figNum,solution)
    %% ########################################################################
    %  FUNCTION TEST_INITIAL_CONDITIONS
    %  
    %  Version: 1.0
    %
    %  INPUT:
    %
    %  OUTPUT:
    %
    %  HISTORY:
    %      v1.0: First runable state. (2014/07/21)
    %
    %% ########################################################################
    
    %% INITILIZATION
    %--------------------------------------------------------------------------
    format shortg;
    %format compact;
    set(0,'DefaultFigureWindowStyle','docked')

    % Change the current folder to the folder of this m-file.
    if(~isdeployed)
        cd(fileparts(which(mfilename)));
    end
    
    %%
    
    
    %% Too make function definitions easier we define a bunch of global variables

    global I_Sk l_Sk p_Sk w_Sk
    global fit_I_Sk fit_l_Sk fit_p_Sk fit_w_Sk
    global Int_F
    global p order

    %%


    %% Code

    % This is a try to see the effect of different starting values
    ord1 = floor(log10(abs(poly(1))));
    ord2 = floor(log10(abs(poly(2))));
    X = 10 .^ (1:ord1);
    Y = 10 .^ (1:ord2);
    Z = zeros(length(Y),length(X));
    max_x = 0;
    max_y = 0;
    max_sol = [];
    max_val = 1;
    min_func = @(x)observe_compensation_min_Func(x,omega,intensity,phase,lambda',fourierlimit);
    options = optimset('MaxIter', 1000,'MaxFunEvals',1e5);
    poly_safe = poly;
    for x=1:length(X)
        for y=1:length(Y)
            p = poly(1:end-2);
            p(1) = X(x);
            p(2) = Y(y);
            [sol2,val2] = fminsearch(min_func,p,options);
            if val2 < max_val
                max_val = val2;
                max_x = X(x);
                max_y = Y(y);
                max_sol = sol2;
            end
            Z(y,x) = val2;
        end
    end
    
    % Ok we found an optimum. Now adjust the linear and constant term.
    % Don't forget to append back the last two constants at the end
    max_sol = [max_sol poly(end-1:end)];
    
    % Now we have to adjust the linear and constant term so that we can see
    % a comparison with the other phases
    min_func = @(x)sum((polyval(poly,omega)-polyval(max_sol,omega)-polyval(x,omega)).^2);
    x0 = [0,0];
    [s,v] = fminsearch(min_func,x0);
    max_sol(end-1:end) = max_sol(end-1:end) + s;
    
    
    fprintf('The best result: %2.2f is won for (%d,%d)\n',max_val,max_x,max_y)
    fprintf('Final constants:\n')
    disp(max_sol)
    figure(figNum)
    disp(Z)
    surf(X,Y,Z)
    figure(figNum+1)
    hold on
    plot(omega,phase,'k')
    P = polyval(poly,omega);
    plot(omega,P,'b')
    P = polyval(solution,omega);
    plot(omega,P,'r')
    P = polyval(max_sol,omega);
    plot(omega,P,'g')
    hold off
    legend('original','least square','least square guess','range guess')
    disp(poly)
    disp(solution)
    disp(max_sol)
end