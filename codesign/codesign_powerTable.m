% Tabulates performance for complex-conjugate (CC - perfect impedance
% matching) and proportional-integral (PI) controllers designed based on
% either the hydro-mechanical system (e.g., "CC on mech") or the electrical
% system (e.g., "CC on elec") in a single sea state.

% Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC
% (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
% Government retains certain rights in this software.
%
% This file is part of fbWecCntrl.
%
%     fbWecCntrl is free software: you can redistribute it and/or modify it
%     under the terms of the GNU General Public License as published by the
%     Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
%
%     fbWecCntrl is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with fbWecCntrl.  If not, see <https://www.gnu.org/licenses/>.
% ---------------------------------------------------------------------

%% Load WEC device data

cf = 60;
mf = load('waveBot_heaveModel.mat');
Zi = mf.Zi_frf(cf:end,1);
Hex = mf.H_frf(cf:end,1)*1e1;
f = mf.f(cf:end,1);
w = 2*pi*f;
dw = w(2)-w(1);

Zpto = PTO_Impedance([1, 0, 0, 0, sqrt(2/3), 1e-3, 0],w); % [N, Id, Bd, Kd, Kt, Rw, Lw]

%% Define sea state and excitation

Hs = 0.125;
gamma = 3.3;
Tp = [1.58, 2.5, 3.5];
for ii = 1:length(Tp)
    S = jonswap(w, [Hs, Tp(ii), gamma]);    % Wave energy density spectrum
    A = sqrt(2*dw*S.S(:));              % wave amplitude spectrum
    Fe(:,ii) = A .* Hex(:);
    Fe_name{ii} = sprintf('Tp%.2f',Tp(ii));
end

%% Design controllers

for ii = 1:size(Fe,2)
    fprintf('\n\n%s\n------------------------\n',Fe_name{ii})
    [T{ii},Pmax(:,ii),mPmech(:,:,ii),mPelec(:,:,ii)] = myfunc(Fe(:,ii),Zi,Zpto,w);
    vn = T{ii}.Properties.VariableNames;
    for jj = 1:length(vn)
        vn{jj} = sprintf('%s_%.2f',vn{jj},Tp(ii));
    end
    T{ii}.Properties.VariableNames = vn;
end

%% Plotting

nf = max(Pmax);
fa = 0.125;

fig = figure;
fig.Position = fig.Position .* [1 1 1 0.75];

t = tiledlayout(2,1);
t.TileSpacing = 'compact';

ax(1) = nexttile;
hold on
grid on

for ii = 1:size(Pmax,2)
    plt.pm(ii) = area(ax(1), f, Pmax(:,ii) / nf(ii), 'facealpha', fa);
end

ax(2) = nexttile;
hold on
grid on

for ii = 1:size(Pmax,2)
    plt.pi(ii) = plot(ax(2), f, mPelec(:,4,ii) ./ mPelec(:,3,ii),...
        'LineWidth',1.5);
    lgs{ii} = sprintf('$T_p = %.2fs$',Tp(ii));
end

ylabel(ax(1),'Norm. $P^{max}_m$ [ ]','interpreter','latex')

legend(ax(1),lgs,'interpreter','latex')
ylabel(ax(2),'$P^\prime_L$ [ ]','interpreter','latex')
xlabel(ax(2),'Frequency [Hz]','interpreter','latex')
% set(ax(2),'yscale','log')

linkaxes(ax,'x')
set(ax(1),'xticklabel',[])
xlim([0.2, 1])

exportgraphics(fig,'codesign_powerRatio.pdf','ContentType','vector')

%% Inline function

function [T,Pmax,mPmech,mPelec] = myfunc(Fe,Zi,Zpto,w)
    
    optimOpts = optimoptions('fminunc',...
        'MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');
    
    Pmax = abs(Fe).^2 ./ (8*real(Zi));
    
    %---------------------------------
    wc(1).leg = 'CC on mech';
    wc(1).ZL = Zi2ZL(Zpto, conj(Zi));
    
    %---------------------------------
    wc(2).leg = 'PI on mech';
    wc(2).cinfo.type = 'PI';
    wc(2).cinfo.w = w;
    wc(2).cinfo.x0 = ones(1,2)*0.1;
    wc(2).objfun = @(x) Pmech( Zi2ZL(Zpto,fbc(x,wc(2).cinfo)),...
        Zpto,...
        Zi,Fe );
    [wc(2).y, wc(2).fval] = fminunc(wc(2).objfun, wc(2).cinfo.x0, optimOpts);
    wc(2).ZL = Zi2ZL(Zpto,fbc(wc(2).y, wc(2).cinfo));
    
    %---------------------------------
    wc(3).leg = 'CC on elec';
    wc(3).ZL = conj( squeeze(Zpto(2,2,:)) ...
        - squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:)) ...
        ./ (squeeze(Zpto(1,1,:)) + Zi) );
    
    %---------------------------------
    wc(4).leg = 'PI on elec';
    wc(4).cinfo.type = 'PI';
    
    wc(4).cinfo.w = w;
    wc(4).cinfo.x0 = ones(1,2);
    wc(4).objfun = @(x) Pelec( Zi2ZL(Zpto,fbc(x,wc(4).cinfo)),...
        Zpto,...
        Zi,Fe );
    [wc(4).y, wc(4).fval] = fminunc(wc(4).objfun, wc(4).cinfo.x0, optimOpts);
    wc(4).ZL = Zi2ZL(Zpto,fbc(wc(4).y, wc(4).cinfo));
    
    
    %% Calculate power and efficiency
    
    for ii = 1:length(wc)
        
        % evaluate mech performance
        [Pmech_tot(ii), mPmech(:,ii)] = Pmech(wc(ii).ZL, Zpto, Zi, Fe);
%         assert(-1*Pmech_tot(ii) <= sum(Pmax),...
%             sprintf('''%s'' making more mechanical power than theoretical limit',wc(1).leg))
        
        % evaluate mech performance
        [Pelec_tot(ii), mPelec(:,ii)] = Pelec(wc(ii).ZL, Zpto, Zi, Fe);
%         assert(-1*Pelec_tot(ii) <= sum(Pmax),...
%             sprintf('''%s'' making more mechanical power than theoretical limit',wc(1).leg))
        
        legCel{ii} = wc(ii).leg;
    end
    
    eta_mech = Pmech_tot./(-1 * sum(Pmax));
    eta_elec = Pelec_tot./(-1 * sum(Pmax));
    
    %% Tabulate results
    
    T = table(-1*Pmech_tot(:)/1e3,eta_mech(:),-1*Pelec_tot(:)/1e3,eta_elec(:),...
        'VariableNames',{'MechPow_kW','MechEfficiency','ElecPow_kW','ElecEffciency'},...
        'RowNames',legCel);
    disp(T)
    
end
