function [type,B,C]=get_sellmeier_coeff(material,varargin)
    if (nargin==1)
       polarization='o'; 
    elseif (nargin==2)
        polarization=varargin{1};
    elseif (nargin>2)
            error('Wrong number of input arguments in "get_sellmeier_coeff"');
    end
    
    
    switch material
        case 'MgF2'
            if ( polarization=='e')
                B=[4.13440230e-1 5.04974990e-1 2.49048620];
                C=[1.35737865e-3 8.23767167e-3 5.65107755e2];
            else
                B=[4.87551080e-1 3.98750310e-1 2.31203530];
                C=[1.88217800e-3 8.95188847e-3 5.66135591e2];
            end
            type=1;
        case 'Sapphire'
            if ( polarization=='e')
                B=[1.50397590e-0 5.50691410e-1 6.59273790];
                C=[5.48041129e-3 1.47994281e-2 4.02895140e2];
            else
                B=[1.43134930e-0 6.50547130e-1 5.34140210];
                C=[5.27992610e-3 1.42382647e-2 3.25017834e2];
            end
            type=1;
        case 'CaF2'
            B=[5.67588800e-1 4.71091400e-1 3.84847230];
            C=[2.52642999e-3 1.00783328e-2 1.20055597e3];
            type=1;
        case 'FS'
            B=[6.96166300e-1 4.07942600e-1 8.97479400e-1];
            C=[4.67914826e-3 1.35120631e-2 9.79340025e1];
            type=1;
        case 'Schott BK7'
            B=[1.03961212e-0 2.31792344e-1 1.01046945e-0];
            C=[6.00069867e-3 2.00179144e-2 1.03560653e2];
            type=1;
        case 'YAG'
            B=[2.2779e-0 0 0];
            C=[1.142e-2 0 0];
            type=1;
        case 'BBO'
            if ( polarization=='e')
                B=[2.3753 0.01224 -0.01516];
                C=[0.01667];   
            else
                B=[2.7359 0.01878 -0.01354];
                C=[0.01822];
            end
            type=2;
%         case 'Schott N-BK7'
%         case 'Schott F2'
        case 'Schott N-LASF46B'
            B=[2.17988922 0.306495184 1.56882437];
            C=[0.0125805384 0.0567191367 105.316538];
            type=1; 
        case 'Schott N-SF14'
           B=[1.69022361 0.28870052 1.7045187];
           C=[0.0130512113 0.061369188 149.517689];
           type=1; 
        case 'Schott SF6'
           B=[1.72448482 0.390104889 1.04572858];
           C=[0.0134871947 0.0569318095 118.557185];
           type=1; 
%         case 'Schott SF10'
%         case 'Schott N-SF10'
        case 'Schott SF11'
           B=[1.73848403 0.311168974 1.17490871];
           C=[0.0136068604 0.0615960463 121.922711];
           type=1;   
        case 'Schott SF57'
           B=[1.81651371 0.428893641 1.07186278];
           C=[0.0143704198 0.0592801172 121.419942];
           type=1;
%         case 'Schott N-SF11'
        otherwise
            error('glass not in list');
            
    end
            
end