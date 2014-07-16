function [x n0 v_g D]=sellmeier(material,lambda0,la, varargin)
		 % Sellmeier function returns refr. index the GVD
		 % and the group velocity
         % v_g will be returned in m/s
         % D will be returned in s^2/m
         if(nargin == 4)
            polarization=varargin{1};
            low=1;
            high=max(size(la));
         elseif nargin == 5
             polarization=varargin{1};
		    low=varargin{2};
		    high=size(la);
         elseif nargin == 6
             polarization=varargin{1};
		    low=varargin{2};
		    high=varargin{3};
         else
            polarization='o';
            low=1;
            high=max(size(la)); 
         end
            
           [type B C]=get_sellmeier_coeff(material, polarization);
            
           c=299792458;
           lambda=la(low:high).*1e6; %lambda in um
           k=size(la);
           sprintf('Sellmeier B-coefficient: %.8f\n',B);
           sprintf('Sellmeier C-coefficient: %.8f\n',C);     
           sprintf('Carrier wavelength: %.2f nm',lambda0*1e9);
           
           if(type==1)
               n=sqrt(1 + B(1).*lambda.^2./(lambda.^2-C(1)) ...
                   + B(2).*lambda.^2./(lambda.^2-C(2)) ...  
                   + B(3).*lambda.^2./(lambda.^2-C(3))); 
               n0=sqrt(1 + B(1).*(lambda0.*1e6).^2./((lambda0.*1e6).^2-C(1)) ...
                   + B(2).*(lambda0.*1e6).^2./((lambda0.*1e6).^2-C(2)) ...  
                   + B(3).*(lambda0.*1e6).^2./((lambda0.*1e6).^2-C(3)));
           elseif(type==2)
               n=sqrt(B(1) + B(2)./(lambda.^2-C(1)) ...
                   + B(3).*lambda.^2);
               n0=sqrt(B(1) + B(2)./((lambda0.*1e6).^2-C(1)) ...
                   + B(3).*(lambda0.*1e6).^2);
           else
               n=0; n0=0;
           end
           %plot(lambda,n)
           
           
           if (low > 1)
             n1 = ones(1,low-1);
             n=[n1 n];
           end
           if (high < k(2))
             n2 = ones(1,k(2)-high);
             n=[n n2];
           end

           diff_n=diff(n)./diff(la);
           diffdiff_n=diff(diff_n)./(diff(la(1:end-1))); 
           
            
           
           tmp = abs(la-lambda0);
           [C,index] = min(tmp);
           if (index> k(2)-2)
               disp('Index too large')
               index = k(2)-2;
           end
           
           [tmp,indn]   = min(abs(la-lambda0));
           [tmp,indv_g] = min(abs(la-lambda0));
           v            = c.*(n(indn)-lambda0.*diff_n(indv_g)).^(-1);
           
           D_nu = lambda0.^3./c.^2.*diffdiff_n(index);
           sprintf('Sellmeier C-coefficient: %.3f\n fs^2/mm',D_nu*1e27);
           
           n0 = [n0, diff_n(indv_g), diffdiff_n(index)];
           %index=floor(lambda0/(lambda(2)-lambda(1)))+1; %achtung index problem
           %v_g_inv = n(index)/obj.c - 2.*pi.*obj.c./lambda0.*diff_n(index);
 
           x={n; diff_n; diffdiff_n};
           v_g = v;
           D=D_nu;
end
