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
%% Clear everything

clc
clear
close all

%% Set up problem

% load device model
mf = matfile(fullfile('data','foswec_model.mat'));
f = mf.f;       % frequency vector [Hz]
Hex = mf.Hex;   % excitation model (in: wave elev., out: flap torque)
Zi = mf.Zi;     % impedance model (in: flap velocity, out: flap torque)

Kt = 0.943;     % motor torque constant
R = 1.082;    % motor electrical winding resistance (0 for mech power)
N = 3.75;       % gear ratio


% replicate for symetric MIMO system
R = eye(size(Hex,1))*R;     
Kt = eye(size(Hex,1))*Kt;
N = eye(size(Hex,1))*N;

% assemble into structure
motspecs = {Kt, R, N};      

% generate sea state
Hm0 = 0.127; 
Tp = 3;
gamma = 1.0;
S = jonswap(2*pi*f,[Hm0, Tp, gamma]);

opts.symFlag = 1;   % force controller to be symmetric
opts.diagFlag = 1;  % allow for cross-coupling control

plotFlag = 1;       % enable plotting

%% Run analysis

[powStudy,fh] = runPowStudy(f,Zi,Hex,S,motspecs,plotFlag,opts);
disp(powStudy(1).gainMatrix)    % displays gain matrix as complex matrix
                                % without frequency dependency, such that
                                % k(i,j) = kP(i,j) -1i kI(i,j)

% set xlims
for ii = 1:length(fh)
    figure(fh(ii))
    xlim([0.2,1])
end

return

%% Save figures

export_fig(fullfile('.','foswec_fbController.pdf'),...
    '-p2',fh(3))
export_fig(fullfile('.','foswec_impedanceMatching.pdf'),...
    '-p2','-painters',fh(4))
