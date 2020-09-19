function C = fbc(x,cinfo)
    % fbc   freq. response func. for feedback controller
    
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
        x (1,:) double {mustBeReal, mustBeFinite}
        cinfo (1,1) struct
    end
    
    switch cinfo.type
        case 'PI'
            C = x(1) - 1i*x(2)./cinfo.w(:);
        case 'P'
            C = x(1);
        case 'CC'
            C = conj(cinfo.Zi);
        otherwise
            error('Invalid value for ''cinfo.type''')
    end
end
