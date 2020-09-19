function Z_mat = PTO_Impedance(Gamma, w)
    % PTO_Impedance   freq. response func. power take-off (PTO) impedance
    
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
       Gamma (1,7) double {mustBeReal, mustBeFinite}
       w (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
    end
    
    Nfreq = length(w);
    
    Z_mat = zeros(2,2,Nfreq);
    
    N = Gamma(1);
    Id = Gamma(2);
    Bd = Gamma(3);
    Kd = Gamma(4);
    Kt = Gamma(5);
    Rw = Gamma(6);
    Lw = Gamma(7);
    
    Zd = Bd + 1i*(w(:)*Id - Kd./w(:));
    Zw = Rw + 1i*w(:)*Lw;
    
    Z_mat(1,1,:) = Zd * N^2;
    Z_mat(1,2,:) = -Kt * N * sqrt(3/2);
    Z_mat(2,1,:) = Kt * N * sqrt(3/2);
    Z_mat(2,2,:) = Zw;
    
end
