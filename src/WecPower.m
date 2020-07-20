function [P, P_f, P_ub, Pub_f]= WecPower(Zi, Fe, C, Kt, R, N)
% WecPower  Calculates power based on impedance, Zi, excitation, Fe, and
% feedback controller, C (along with PTO specifications).
%
% Args:
%   Zi      impedance model (in: velocity, out: torque) with 
%           size(Zi) = [nDof,nDof,nFreq]
%   Fe      excitation spectrum with size(Fe) = [nDof, nFreq]
%   C       control model (in: velocity, out: torque) with 
%           size(C) = [nDof,nDof,nFreq]
%   Kt      torque constant with size(Kt) = [nDof, nDof] 
%           and isdiag(Kt) = 1
%   R       motor winding resistance with size(R) = [nDof, nDof] 
%           and isdiag(R) = 1
%   N       gear ratio (in: torque at motor, out: torque at body) with
%           size(N) = [nDof, nDof] and isdiag(N) = 1

% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------

Ke = Kt*2/3;                    % 3-phase PMS motor elect. const.

nFreq = size(Zi,3);

% preallocate arrays 
Pfh = zeros(nFreq,1);           % complex power
Pub_f = zeros(nFreq,1);         % upper bound ('complex conjugate control')
Omega = zeros(length(R),nFreq);

for ii = 1:nFreq
    
    Omega(:,ii) = ( Zi(:,:,ii) - C(:,:,ii) ) \ Fe(:,ii);
    
    Pfh(ii) = ( (Ke*N + (R/(N*Kt)) * C(:,:,ii)) * Omega(:,ii) )'...
              * ( ((N*Kt) \ C(:,:,ii)*Omega(:,ii)) );
          
    Pub_f(ii) = -1 * real(1/8 * (Fe(:,ii)' / real(Zi(:,:,ii))) * Fe(:,ii));
end

P_f = real(Pfh) * 3/4;
P_ub = sum(Pub_f);
P = sum(P_f);

end