function [P, P_f, Z_L, den_P] = co_optimize_oneDOF_PTO_power(x, PTO_param_opt_mask, PTO_params, w, Zi, Fe)

n_PTO_var = sum(PTO_param_opt_mask);
if n_PTO_var > 0
    PTO_params(PTO_param_opt_mask) = x(1:n_PTO_var);
end

Z_L = Load_impedance(w, x(n_PTO_var+1:end));

Z_PTO = PTO_Impedance(w, PTO_params);

den_P = ((squeeze(Z_PTO(2,2,:)) + Z_L) .*...
    ( squeeze(Z_PTO(1,1,:)) + Zi )) - squeeze(Z_PTO(2,1,:).*Z_PTO(1,2,:) ) ;
P_f = 0.5 * real(Z_L) .* (abs((Fe .* squeeze(Z_PTO(2,1,:)) ) ).^2 ./ abs(den_P).^2);
P = -sum(P_f);



% N = PTO_params(1);
% Id = PTO_params(2);
% Bd = PTO_params(3);
% Kd = PTO_params(4);
% Kt = PTO_params(5);
% Rw = PTO_params(6);
% Lw = PTO_params(7);
% 
% 
% Zd = Bd + 1i*(w(:)*Id - Kd./w(:));
% Zw = Rw + 1i*w(:)*Lw;



% den_P = ( (((Zd+Zi/N^2)./(Kt*sqrt(3/2) )) .* (1 + Zw./Z_L)) -(sqrt(3/2)*Kt ./ Z_L ) );
% P_f = abs( (Fe/N) ./ den_P ).^2 ./ (2*real(Z_L));
% 
% den_P = (Zw+Z_L).*(N^2*Zd+Zi) - (Kt*N)^2 *2/3 ;
% P_f = (abs(( Kt*N*Z_L.*Fe )./( den_P )).^2) ./ (2*real(Z_L));
% 
% P = -sum(P_f);
end