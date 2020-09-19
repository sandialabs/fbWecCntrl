function [P, P_f, Z_L, Z_PTO, C] = co_optimize_oneDOF_PTO_power(x, PTO_param_opt_mask, PTO_params, w, Zi, Fe, Cfunc)
    % co_optimize_oneDOF_PTO_power co-optimization of PTO and control

    % Copyright 2020 National Technology & Engineering Solutions of Sandia,
    % LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the
    % U.S. Government retains certain rights in this software.
    %
    % This file is part of fbWecCntrl.
    %
    %     fbWecCntrl is free software: you can redistribute it and/or
    %     modify it under the terms of the GNU General Public License as
    %     published by the Free Software Foundation, either version 3 of
    %     the License, or (at your option) any later version.
    %
    %     fbWecCntrl is distributed in the hope that it will be useful, but
    %     WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    %     General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with fbWecCntrl.  If not, see
    %     <https://www.gnu.org/licenses/>.
    % ---------------------------------------------------------------------
    
    arguments
        x (1,:) double {mustBeReal,mustBeFinite}
        PTO_param_opt_mask (1,:) {mustBeNumericOrLogical}
        PTO_params (1,:) double {mustBeReal,mustBeFinite}
        w (:,1) double {mustBeReal,mustBeFinite,mustBePositive}
        Zi (:,1) double {mustBeFinite}
        Fe  (:,1) double {mustBeFinite}
        Cfunc (1,:) string = 'PI'
    end
    
    % PTO impedance -------------------------------------------------------
    n_PTO_var = sum(PTO_param_opt_mask);
    if n_PTO_var > 0
        PTO_params(PTO_param_opt_mask) = x(1:n_PTO_var);
    end
    
    Z_PTO = PTO_Impedance(PTO_params,w);
    
    % Load impedance ------------------------------------------------------
    switch Cfunc
        case 'P'
            xc = x(end);
        case 'PI'
            xc = x(end-1:end);
        case 'CC'
            xc = [];
        otherwise
            error('Invalid value for ''Cfunc''')
    end
    
    cinfo.w = w;
    cinfo.Zi = Zi;
    cinfo.type = Cfunc;
    C = fbc(xc,cinfo);

    Z_L = Zi2ZL(Z_PTO, C);
    
    % power at load -------------------------------------------------------
    [P, P_f] = Pelec(Z_L, Z_PTO, Zi, Fe);
    
end
