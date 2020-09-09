function [Pmech_tot, Pmech] = Pmech(ZL, Zpto, Zi, Fe)
    
    Zin = squeeze(Zpto(1,1,:)) ...
        - ((squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:))) ...
        ./ (squeeze(Zpto(2,2,:)) + ZL));
    
    Pmech = 1/2 * abs( Fe ./ (Zi + Zin) ).^2 .* real(Zin);
    
    Pmech_tot = -1 * sum(Pmech);
end
