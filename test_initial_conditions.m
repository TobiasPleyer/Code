function test_initial_conditions(solution)
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

    global w_Sk p_Sk
    global p figNum

    %%


    %% Code

    % This is a try to see the effect of different starting values
    ord1 = floor(log10(abs(p(1))));
    ord2 = floor(log10(abs(p(2))));
    X = floor(p(1)/(10^(ord1-1)))-5:0.4:floor(p(1)/(10^(ord1-1)))+5;
    X=X.*(10^(ord1-1));
    Y = floor(p(2)/(10^(ord2-1)))-5:0.4:floor(p(2)/(10^(ord2-1)))+5;
    Y=Y.*(10^(ord2-1));
    Z = zeros(length(Y),length(X));
    max_x = 0;
    max_y = 0;
    max_sol = [];
    max_val = 1;
    min_func = @(x)observe_compensation_min_Func(x);
    options = optimset('MaxIter', 1000,'MaxFunEvals',1e5);
    for x=1:length(X)
        for y=1:length(Y)
            p_new = p(1:end-2);
            p_new(1) = X(x);
            p_new(2) = Y(y);
            [sol2,val2] = fminsearch(min_func,p_new,options);
            if val2 < max_val
                max_val = val2;
                max_x = X(x);
                max_y = Y(y);
                max_sol = sol2;
            end
            Z(y,x) = val2;
        end
    end
    % Find out what the best five values are
    [sortedValues,sortIndex] = sort(Z(:),'descend');
    maxIndex = sortIndex(end-4:end);
    fprintf('Best values:\n')
    disp(Z(maxIndex))
    
    % Ok we found an optimum. Now adjust the linear and constant term.
    % Don't forget to append back the last two constants at the end
    max_sol = [max_sol p(end-1:end)];
    
    % Now we have to adjust the linear and constant term so that we can see
    % a comparison with the other phases
    min_func = @(x)sum((polyval(p,w_Sk)-polyval(max_sol,w_Sk)-polyval(x,w_Sk)).^2);
    x0 = [0,0];
    [s,v] = fminsearch(min_func,x0);
    max_sol(end-1:end) = max_sol(end-1:end) + s;
    
    
    fprintf('The best result: %2.2f is won for (%d,%d)\n',max_val,max_x,max_y)
    fprintf('Final constants:\n')
    disp(max_sol)
    figure(figNum)
    figNum = figNum + 1;
    disp(Z)
    surf(X,Y,Z)
    title('3D visualization of the optimization outcome for different initial conditions around our least square guess')
    xlabel('Highest order polynomial constant')
    ylabel('Second highest order polynomial constant')
    zlabel('Achieved optimum')
    figure(figNum)
    figNum = figNum + 1;
    hold on
    plot(w_Sk,p_Sk,'k')
    P = polyval(p,w_Sk);
    plot(w_Sk,P,'b')
    P = polyval(solution,w_Sk);
    plot(w_Sk,P,'r')
    P = polyval(max_sol,w_Sk);
    plot(w_Sk,P,'g')
    hold off
    legend('original','least square','least square guess','range guess')
%     disp(p)
%     disp(solution)
%     disp(max_sol)
%     fprintf('Difference of initial value to least square:\n')
%     fprintf('\tHighest order:\n')
%     disp(X-ones(size(X)).*p(1))
%     fprintf('\tSecond highest order:\n')
%     disp(Y-ones(size(Y)).*p(2))
end