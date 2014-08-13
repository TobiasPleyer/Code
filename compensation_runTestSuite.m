function compensation_runTestSuite()
    fprintf('---------RUNNING TEST SUITE---------\n')
    fprintf('l0: 1030 nm\n')
    fprintf('w : 1.82:0.0001:1.85 * 1e15 s^-1\n')
    fprintf('I : exp(-(w-w0).^2/1e-5)\n')
    fprintf('\n')
    fprintf('\n')
    c  = 299792458;
    l0 = 1030;
    w0 = 1e-15*  2*pi*c / (l0*1e-9);
    w  = 1.82:0.0001:1.84;
    l  = 1e9*    2*pi*c ./ (w*1e15);
    p  = 10000*(w-w0).^2; 
    I  = exp(-(w-w0).^2/1e-5);
    [Int_F,t_F,Ek_F] = compensation_calcFourierlimit(I,l);
end