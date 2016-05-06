function RES = resFUN_vdp(x)
%% <Importing global variables>
    global freq
    global func
%% <FFT for given x to solve for derivatives xddot and xdot>
    X = fft(x);
    xdot = ifft(1i*freq.*X);
    xddot = ifft(-freq.^2.*X);
    
%% <Identifying residual function>
    res = xddot+xdot.*(x.^2-1)+x-func;
    RES = sum(abs(res.^2));

    