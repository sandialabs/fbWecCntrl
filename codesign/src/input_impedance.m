function Zin = input_impedance(Zpto, ZL)
    % input_impedance   input impedance freq. response func.
    
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
        Zpto (2,2,:) double {mustBeFinite}
        ZL (:,1) double {mustBeFinite}
    end
    
    Zin = squeeze(Zpto(1,1,:)) ...
        - squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:)) ...
        ./ (squeeze(Zpto(2,2,:)) + ZL);
end