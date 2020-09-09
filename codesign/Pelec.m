function [Pelec_tot, Pelec] = Pelec(ZL, Zpto, Zi, Fe)
    
    Pelec = 0.5 *real(ZL) .* abs( ( squeeze(Zpto(2,1,:)) .* Fe(:)) ./ ...
        ( (( squeeze(Zpto(1,1,:)) + Zi ) .* ( squeeze(Zpto(2,2,:)) + ZL )) ...
        - squeeze(Zpto(1,2,:) .* Zpto(2,1,:)) ) ).^2;
    
    Pelec_tot = -1 * sum(Pelec);
    
end