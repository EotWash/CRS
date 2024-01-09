function noise=CRSTheoryNoiseModel(f)

%     Q=18.2;
    Q=1000;
    f0=16.7e-3;
    r=15e-2; 
    kb=1.380e-23; 
    T=293;    
    I=0.0941;
    
    Tf=1./abs(f.^2./(f.^2-i*(f0^2/Q)-f0^2));

    readout=Tf.*2e-13/r;
    thermal = Tf.*sqrt(4*kb*T*f0^2./(I*2*pi*f.*Q.*((f.^2-f0^2).^2+f0^4/Q^2)));

    noise=sqrt(readout.^2+thermal.^2);
    
end