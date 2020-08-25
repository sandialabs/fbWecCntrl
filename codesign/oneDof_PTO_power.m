function [P, P_f] = oneDof_PTO_power(Z_L, Z_pto, Zi, Fe)
    
    den_P = ((squeeze(Z_pto(2,2,:)) + Z_L) .* ...
        ( squeeze(Z_pto(1,1,:)) + Zi )) - squeeze(Z_pto(2,1,:).*Z_pto(1,2,:) ) ;
    P_f = 0.5 * real(Z_L) .* abs((Fe .* squeeze(Z_pto(2,1,:)) ) ./ den_P ).^2 ;
    P = -sum(P_f);
    
end