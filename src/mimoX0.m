function [x0,nVars] = mimoX0(nGains,nDof,symFlag,diagFlag)
% Finds the number of gains for a given MIMO feedback controller and
% creates a column vector of zeros for an initial guess.
%
% Args:
%   Gains       number of gains (P: 1, PI: 2)
%   nDof        number of degrees of freedom
%   symFlag     set to one for symmetric 
%               (M == M.' AND rot90(M) == rot90(M.')
%   diagFlag    set to one for diagonal (isdiag(M) == 1)

% Copyright 2020 National Technology & Engineering Solutions of Sandia,
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the
% U.S. Government retains certain rights in this software.
%
% This file is part of fbWecCntrl.
%
%     fbWecCntrl is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     fbWecCntrl is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with fbWecCntrl.  If not, see <https://www.gnu.org/licenses/>.
%
% -------------------------------------------------------------------------

if nDof == 1                    % SISO
    nVars = nGains;
else
    if diagFlag
        if symFlag
            nVars = (floor(nDof/2) + mod(nDof,2)) * nGains;
        else
            nVars = nDof * nGains;
        end
    else
        if symFlag
            nVars = nGains * sum(sum(triu(ones(nDof,nDof),1) ...
                + diag([ones(nDof/2,1);zeros(nDof/2,1)])));
        else
            nVars = nGains * nDof^2;
        end
    end
end

x0 = zeros(nVars,1);

end
