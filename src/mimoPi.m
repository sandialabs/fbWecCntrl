function [C] = mimoPi(x,nGains,nDof,w,symFlag,diagFlag)
% Produces the transfer function for a P or PI controller
%
% Args.
%   x       gains given as a column vector; if P controller, x = [kP]; if P
%           controller, x = [kP;kI];
%   nGains  number of gains (P: 1, PI: 2)
%   nDof    numer of degrees of freedom in the MIMO system (assumed square)
%   w       frequency vector [rad/s]
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

arguments
    x (:,1) {mustBeReal}
    nGains (1,1) {mustBeReal,mustBePositive}
    nDof (1,1) {mustBeReal,mustBePositive}
    w (:,1) {mustBeReal,mustBePositive}
    symFlag (1,1) {mustBeNumericOrLogical}
    diagFlag (1,1) {mustBeNumericOrLogical}
end

% prepare inputs
w = w(:);
nFreq = length(w);
n = length(x);

% if diagFlag && symFlag
%     nVar = floor(nDof/2) + mod(nDof,2);
% end

if nGains == 2
    kP = x(1:n/2);
    kI = x(n/2+1:n);
else
    kP = x;
    kI = 0 * x;
end

if symFlag && nDof > 1
    kP = [kP;flipud(kP)];
    kI = [kI;flipud(kI)];
end

if diagFlag
    CP = repmat(diag(kP),[1,1,nFreq]);
    tmp = diag(kI);
    CI = -1i * reshape(tmp(:)./w',[nDof,nDof,nFreq]);
else
    CP = repmat(reshape(kP,[nDof,nDof]),[1,1,nFreq]);
    CI = -1i * reshape(kI./w',[nDof,nDof,nFreq]);
end

C = CP + CI;

end
