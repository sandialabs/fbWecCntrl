function [P, P_f, den_P] = oneDof_PTO_power(Z_L, PTO_params, w, Zi, Fe)
    
    
    Z_PTO = PTO_Impedance(w, PTO_params);
    
    den_P = ((squeeze(Z_PTO(2,2,:)) + Z_L) .*...
( squeeze(Z_PTO(1,1,:)) + Zi )) - squeeze(Z_PTO(2,1,:).*Z_PTO(1,2,:) ) ;
    P_f = 0.5 * real(Z_L) .* (abs((Fe .* squeeze(Z_PTO(2,1,:)) ) ).^2 ./ abs(den_P).^2);
    P = -sum(P_f);
    
end