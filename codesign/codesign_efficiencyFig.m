% Finds the frequency dependent efficiency for complex-conjugate (CC -
% perfect impedance matching) and proportional-integral (PI) controllers
% designed based on either the hydro-mechanical system (e.g., "CC on mech")
% or the electrical system (e.g., "CC on elec").

clc
clear
close all

%% Load WEC device data

cf = 60;
mf = load('waveBot_heaveModel.mat');
Zi = mf.Zi_frf(cf:end,1);
Hex = mf.H_frf(cf:end,1)*1e1;
f = mf.f(cf:end,1);
w = 2*pi*f;
dw = w(2)-w(1);

Zpto = PTO_Impedance(w,[1, 0, 0, 0, sqrt(2/3), 1e-3, 0]); % [N, Id, Bd, Kd, Kt, Rw, Lw]

%% Wave excitation

Fe = sqrt(8*real(Zi))*1; % const. power excitation

%% Design controllers

% CC on hydro
legCel{1} = 'CC on hydro';
C{1} = conj(Zi);
ZL{1} = Zi2ZL(Zpto,C{1});

% CC on full sys. (eq. 22), output matching condition only
legCel{2} = 'CC on full sys.';
C{2} = conj( squeeze(Zpto(2,2,:)) ...
    - squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:)) ...
    ./ (squeeze(Zpto(1,1,:)) + Zi) );
ZL{2} = C{2};

%% Calculate power and efficiency

Pmax = abs(Fe).^2 ./ (8*real(Zi));

for ii = 1:length(C)
    Zin{ii} = input_impedance(Zpto,ZL{ii});
    Pmech(:,ii) = oneDof_mech_power(Zi, Zin{ii}, Fe);
    [~,Pelec(:,ii)] = oneDof_PTO_power(ZL{ii},Zpto,Zi,Fe);
end

%% Plot results

fig = figure;
fig.Position = fig.Position .* [1 1 1 0.5];
hold on
grid on

for ii = 1:length(C)
    plt(1,ii) = plot(f,Pmech(:,ii),'--','LineWidth',1.5);
end
ax = gca; ax.ColorOrderIndex = 1;
for ii = 1:length(C)
    plt(2,ii) = plot(f,Pelec(:,ii),'-','LineWidth',1.5);
end

l1 = legend([plt(2,:)],[legCel(:)']);
set(l1,'location','southwest')
ylim([-5,1])
xlim([0.2, 1])

ylabel('Efficiency [ ]')
xlabel('Frequency [Hz]')

% exportgraphics(gcf,'codesign_freqPowerComp.pdf','ContentType','vector')

%% Copyright

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
