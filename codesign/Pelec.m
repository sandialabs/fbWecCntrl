function [Pelec_tot, Pelec] = Pelec(ZL, Zpto, Zi, Fe)
    % Pelec   electrical power
    
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
        ZL (:,1) double {mustBeFinite}
        Zpto (2,2,:) double {mustBeFinite}
        Zi (:,1) double {mustBeFinite}
        Fe (:,1) double {mustBeFinite}
    end
    
    Pelec = 0.5 *real(ZL) .* abs( ( squeeze(Zpto(2,1,:)) .* Fe(:)) ./ ...
        ( (( squeeze(Zpto(1,1,:)) + Zi ) .* ( squeeze(Zpto(2,2,:)) + ZL )) ...
        - squeeze(Zpto(1,2,:) .* Zpto(2,1,:)) ) ).^2;
    
    Pelec_tot = -1 * sum(Pelec);
    
end