function [new_phase,grating_phase,GDD_out,TOD_out,varargout] = compensation_applyGrating(omega,phase,varargin)
    % The formulas used below are taken from the script of Prof. Dr. Karsch,
    % pages 31-34 and apply to a grating compressor.

    strUsage = 'Usage:\napplyGrating(omega,phase)\nor\napplyGrating(omega,phase,G,AOI,N)\n';
    strExplanation = 'Where\nomega: Angular frequency vector\nphase: Phase in angular frequency domain\nG: Distance of the gratings [m] (perpendicular)\nAOI: Angle of incidence [rad]\nN: Line number of the grating (d^-1) [#/mm]\n';
    
    
    % define defaults at the beginning of the code so that you do not need to
    % scroll way down in case you want to change something or if the help is
    % incomplete
    options = struct('G',0.15 ...
                    ,'AOI',2*pi * 45/360 ...
                    ,'N',300e3 ...
                    ,'lambda0', 1.03e-6 ...
                    ,'optimize',false ...
                    );
    
    %# read the acceptable names
    optionNames = fieldnames(options);
    
    if nargin < 2
        fprintf(strUsage)
        fprintf(strExplanation)
    end
    if nargin > 2
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
           error('EXAMPLE needs propertyName/propertyValue pairs')
        end
        for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
            %inpName = lower(pair{1}); %# make case insensitive
            inpName = pair{1};

            if any(strmatch(inpName,optionNames))
                %# overwrite options. If you want you can test for the right class here
                %# Also, if you find out that there is an option you keep getting wrong,
                %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
                options.(inpName) = pair{2};
            else
                error('%s is not a recognized parameter name',inpName)
            end
        end
    end
    
    G = options.G;
    AOI = options.AOI;
    N = options.N;
    d = 1/N;
    lambda0 = options.lambda0;
    
    c = 299792458;
%     lambda = 2*pi*c ./ omega;
    
%     beta = asin(N.*lambda ./ sin(AOI));
%     theta = AOI - beta;
%     
%     full_grating_phase = 2*G .* (omega./c.*((1+cos(theta)./cos(beta))) - tan(beta).*(2*pi)./d);
    
    % Now find the GDD and TOD of our grating
    omega0 = 2*pi*c / lambda0;
    beta0 = asin(N*lambda0-sin(AOI));
    if (nargout == 5)
        if options.optimize
            phase_appr = polyfit(omega,phase,4);
            GDD_orig = polyval(polyder(polyder(phase_appr)),omega0);
            G_optim = fsolve(@(G)1e30.*(-2*(G./cos(beta0))./c.*lambda0./(2*pi*c) .* (lambda0.*N./cos(beta0)).^2+GDD_orig),G);
            GD = G_optim./cos(beta0);
            GDD = -2*GD./c.*lambda0./(2*pi*c) .* (lambda0.*N./cos(beta0)).^2;
            TOD = -3.*GDD.*lambda0./(2*pi*c) .* (1 + lambda0.*N.*sin(beta0)./(cos(beta0)).^2);
            GDD_out = GDD / 1e-30;
            TOD_out = TOD / 1e-45;
            varargout{1} = G_optim;

            grating_phase = TOD./6.*(omega-omega0).^3 + GDD./2.*(omega-omega0).^2;
        else
            GD = G./cos(beta0);
            GDD = -2*GD./c.*lambda0./(2*pi*c) .* (lambda0.*N./cos(beta0)).^2;
            TOD = -3.*GDD.*lambda0./(2*pi*c) .* (1 + lambda0.*N.*sin(beta0)./(cos(beta0)).^2);
            GDD_out = GDD / 1e-30;
            TOD_out = TOD / 1e-45;

            grating_phase = TOD./6.*(omega-omega0).^3 + GDD./2.*(omega-omega0).^2;
            varargout{1} = -1;
        end
    else
        GD = G./cos(beta0);
        GDD = -2*GD./c.*lambda0./(2*pi*c) .* (lambda0.*N./cos(beta0)).^2;
        TOD = -3.*GDD.*lambda0./(2*pi*c) .* (1 + lambda0.*N.*sin(beta0)./(cos(beta0)).^2);
        GDD_out = GDD / 1e-30;
        TOD_out = TOD / 1e-45;

        grating_phase = TOD./6.*(omega-omega0).^3 + GDD./2.*(omega-omega0).^2;
    end
    new_phase = phase + grating_phase;
end