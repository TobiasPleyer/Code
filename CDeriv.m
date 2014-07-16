%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Written by Moritz Ueffing 2014
%       Version 1.0
%
%       calculates nth central derivative evaluating 
%       the function CDervative n times
%       
%       INPUT:  x:      x-vector
%               y:      y-vector
%               n:      Derivative (supported only up to 4th)
%                       
%                   x,y don't need equal spacing as interp1 is 
%                   used in between
%
%       OUTPUT: x_n:    x-vector
%               y_n:    y-vector of nth derivative
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_n, y_n] = CDeriv(x, y, n)
   N = length(y);
   
   xi = linspace(x(1), x(end), N);
   yi = interp1(x, y, xi, 'spline')';
   
   x_p = xi;
   y_p = yi;
   
   if n<=4
    [x_p, y_p] = CDerivative(x_p,y_p,n);
   else
       error('ERROR: Derivative >4 not supported by CDeriv!!');
   end
  
   
   x_n = x_p;
   y_n = y_p';         
end

function [x, f] = CDerivative(x,y,n)
    c = GetCoefficients_CEN(n);
    
    h = abs(x(2)-x(1));
    N = length(y);
    i=1:N;
   
    %construct matrix here
    A = sparse(N,N,0);
    k=length(c)-1;
    for j = 0:k
        shift = j-k/2;
        if (shift <= 0)
            A = A + sparse(i(1-shift:end), i(1:end+shift),c(j+1),N,N);
        else
            A = A + sparse(i(1:end-shift),i(shift+1:end),c(j+1),N,N);
        end
    end
    
    %Initialize first points with FWRD Difference
    c_fwrd = GetCoefficients_FWRD(n);
    length_cfrwd = length(c_fwrd)-1;
    for j = 1:k/2
        A(j,:) = 0;
        A(j,j:j+length_cfrwd) = c_fwrd;
    end
    
    %Initialize last points with BWRD Difference
    c_bwrd = GetCoefficients_BWRD(n);
    length_cbrwd = length(c_bwrd)-1;
    for j = 1:k/2
        A(end-j+1,:) = 0;
        A(end-j+1,end - length_cbrwd-j+1:end-j+1) = c_bwrd;
    end
    
    f = 1/h.^n * A*y;
end

function c = GetCoefficients_CEN(n)
    %O(h^4) coefficients

   switch n
       case 1
            c = [1/12 -2/3 0 2/3 -1/12]; 
       case 2
            c = [-1/12 4/3 -5/2 4/3 -1/12]; 
       case 3
            c = [1/8 -1 13/8 0 -13/8 1 -1/8];
       case 4
            c = [-1/6 2 -13/2 28/3 -13/2 2 -1/6];
       otherwise
           c = [1/12 -2/3 0 2/3 -1/12];
   end
end

function c = GetCoefficients_BWRD(n)
    %O(h^4) coefficients
   c_tmp = GetCoefficients_FWRD(n);
   if (mod(n,2)==1)
       c = -1*c_tmp(end:-1:1);
   else
       c = c_tmp(end:-1:1);
   end
end

function c = GetCoefficients_FWRD(n)
    %O(h^4) coefficients
   switch n
       case 1
            c = [-25/12 4 -3 4/3 -1/4]; 
       case 2
            c = [15/4 -77/6 107/6 -13 61/12 -5/6]; 
       case 3
            c = [-49/8 29 -461/8 62 -307/8 13 -15/8];
       case 4
            c = [28/3 -111/2 142 -1219/6 176 -185/2 82/3 -7/2];
       otherwise
           c = [1/12 -2/3 0 2/3 -1/12];
   end 
end









