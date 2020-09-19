function Zgen = Gen_impedance(genParams, w)
    % Gen_impedance     generator impedance freq. response function
        
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
       genParams (1,4) double {mustBeReal, mustBeFinite}
       w (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
    end
    
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
