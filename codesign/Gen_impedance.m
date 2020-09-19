function Zgen = Gen_impedance(genParams, w)
    
    Ir = genParams(1);      % rotor inertia
    Kt = genParams(2);      % torque coefficient
    Rw = genParams(3);      % winding resistance
    Lw = genParams(4);      % winding inductance
    
    Zw = Rw + 1i*w(:)*Lw;   % winding impedance
    
    nFreq = length(w);
    
    Zgen = zeros(2,2,nFreq);
    
    Zgen(1,1,:) = 1i*w*Ir;
    Zgen(1,2,:) = -sqrt(3/2)*Kt;
    Zgen(2,1,:) = sqrt(3/2)*Kt;
    Zgen(2,2,:) = Zw;
end
