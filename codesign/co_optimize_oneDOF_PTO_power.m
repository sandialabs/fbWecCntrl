function [P, P_f, Z_L, Z_PTO, C] = co_optimize_oneDOF_PTO_power(x, PTO_param_opt_mask, PTO_params, w, Zi, Fe)
    
    % PTO impedance
    % ---------------------------------------------------------------------
    n_PTO_var = sum(PTO_param_opt_mask);
    if n_PTO_var > 0
        PTO_params(PTO_param_opt_mask) = x(1:n_PTO_var);
    end
    
    Z_PTO = PTO_Impedance(w, PTO_params);
    
    % Load impedance
    % ---------------------------------------------------------------------
    
    % PI controller, but can be anything...
    C = x(end-1) - 1i*x(end)./w(:);
    
    Z_L = Load_impedance(Z_PTO, C);
    
    % power at load
    % ---------------------------------------------------------------------
    
    [P, P_f] = oneDof_PTO_power(Z_L, Z_PTO, Zi, Fe);
    
end
