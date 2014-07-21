function test_initial_conditions
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
    %      v1.0: First runable state
    %
    %% ########################################################################
    %% INITILIZATION
    %--------------------------------------------------------------------------
    clc;    % Clear the command window.
    close all;  % Close all figures (except those of imtool.)
%     clear;  % Erase all existing variables. Or clearvars if you want.
    format shortg;
    %format compact;
    set(0,'DefaultFigureWindowStyle','docked')

    % Change the current folder to the folder of this m-file.
    if(~isdeployed)
        cd(fileparts(which(mfilename)));
    end
    %% ########################################################################
    %% Code
    % This is a try to see the effect of different starting values
    fval_history = [];
    iter_history = [];
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
    fprintf('The best result: %2.2f is won for (%d,%d)\n',max_val,max_x,max_y)
    fprintf('Final constants:\n')
    disp(max_sol)
    figure(figNum)
    plot(iter_history,fval_history,'x')
    xlabel('Iteration #')
    ylabel('Returned optimum')
    title(sprintf('Monitoring of the optimization trend curve for order %d.',length(poly)-1))
    figure(22)
    disp(Z)
    surf(X,Y,Z)
    figure(23)
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
end