% Finds the frequency dependent efficiency for complex-conjugate (CC -
% perfect impedance matching) and proportional-integral (PI) controllers
% designed based on either the hydro-mechanical system (e.g., "CC on mech")
% or the electrical system (e.g., "CC on elec").

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

savefigflag = 0;

%% Load WEC device data

cf = 60;
mf = load('waveBot_heaveModel.mat');
Zi = mf.Zi_frf(cf:end,1);
Hex = mf.H_frf(cf:end,1);
f = mf.f(cf:end,1);
w = 2*pi*f;
dw = w(2)-w(1);

Zpto = PTO_Impedance([1, 0, 0, 0, sqrt(2/3), 1e-3, 0],w); % [N, Id, Bd, Kd, Kt, Rw, Lw]

%% Wave excitation

Hs = 0.254;
gamma = 3.3;
Tp = [1.58, 2.5, 3.5];
for ii = 1:length(Tp)
    S = jonswap(w, [Hs, Tp(ii), gamma]);    % Wave energy density spectrum
    A = sqrt(2*dw*S.S(:));              % wave amplitude spectrum
    Fe(:,ii) = A .* Hex(:);
    Fe_name{ii} = sprintf('Tp%.2f',Tp(ii));
end

Fe(:,4) = sqrt(8*real(Zi))*1; % const. power excitation
Fe_name{4} = 'constPowerExc';

Fe(:,5) = ones(size(Fe(:,1)))*mean(Fe(:,1));
Fe_name{5} = 'constFe';

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


for jj = 1:size(Fe,2)
    
    %% Calculate power and efficiency
    
    Pmax = abs(Fe(:,jj)).^2 ./ (8*real(Zi));
    
    for ii = 1:length(C)
        Zin(:,ii) = input_impedance(Zpto,ZL{ii});
        
        V(:,ii) = Fe(:,jj) ./ (Zi + Zin(:,ii));
        Fp(:,ii) = V(:,ii) .* Zin(:,ii);
        X(:,ii) = V(:,ii) ./ (1i * w);
        
        [~,Pmech_f(:,ii)] = Pmech(ZL{ii}, Zpto, Zi, Fe(:,jj));
        
        [~,Pelec_f(:,ii)] = Pelec(ZL{ii}, Zpto, Zi, Fe(:,jj));
    end
    
    %% Plot efficiency
    
    fig(1) = figure('name',sprintf('powerEfficiency_%s',Fe_name{jj}));
    fig(1).Position = fig(1).Position .* [1 1 1 0.5];
    hold on
    grid on
    
    for ii = 1:length(C)
        plt(1,ii) = plot(f,Pmech_f(:,ii) ./ Pmax,'--','LineWidth',1.5);
    end
    ax = gca; ax.ColorOrderIndex = 1;
    for ii = 1:length(C)
        plt(2,ii) = plot(f,Pelec_f(:,ii),'-','LineWidth',1.5);
    end
    
    l1 = legend([plt(2,:)],[legCel(:)']);
    set(l1,'location','southwest')
    ylim([-5,1])
    xlim([0.2, 1])
    
    ylabel('Efficiency [ ]')
    xlabel('Frequency [Hz]')
    
    if savefigflag
        fname = sprintf('codesign_powerEfficiency_%s.pdf',Fe_name{jj});
        exportgraphics(gcf,fname,'ContentType','vector')
    end
    
    %%
    
    fig(2) = figure('name',sprintf('velForceOverdesign_%s',Fe_name{jj}));
    fig(2).Position = fig(2).Position .* [1 1 1 0.75];
    t = tiledlayout(2,1);
    t.TileSpacing = 'compact';
    % t.Padding = 'compact';
    ax(1) = nexttile;
    hold on
    grid on
    set(ax(1),'XTickLabel',[])
    set(ax(1),'YScale','log')
    
    plot(f,abs(X),'LineWidth',1.5);
    legend(legCel)
    ylabel(ax(1),'Motion amp. [m]')
    
    
    ax(2) = nexttile;
    hold on
    grid on
    set(ax(2),'YScale','log')
    plot(f, abs(Fp),'LineWidth',1.5);
    ylabel(ax(2),'PTO force [N]')
    
    xlabel('Frequency [Hz]')
    linkaxes(ax,'x')
    xlim([0.2, 1])
    
    if savefigflag
        fname = sprintf('codesign_velForceOverdesign_%s.pdf',Fe_name{jj});
        exportgraphics(gcf,fname,'ContentType','vector')
    end
    
    %%
    
    fig(3) = figure('name',sprintf('impedances_%s',Fe_name{jj}));
    fig(3).Position = fig(3).Position .* [1 1 1 0.75];
    hold on
    for ii = 1:length(C)
        sys(ii) = frd(Zi + Zin(:,ii),w);
    end
    sys(length(C)+1) = frd(2*real(Zi),w);
    h = bodeplot(sys(1),sys(2),sys(3));
    setoptions(h,'FreqUnits','Hz');
    setoptions(h,'Grid','on');
    setoptions(h,'IOGrouping','all');
    legend(legCel{:},'2*R_i')
    xlim([0.2, 1])
    
    %% Plot power
    
    fig(4) = figure('name',sprintf('power_%s',Fe_name{jj}));
    fig(4).Position = fig(4).Position .* [1 1 1 0.5];
    hold on
    grid on
    
    clear plt
    for ii = 1:length(C)
        plt(1,ii) = plot(f,Pmech_f(:,ii),'--','LineWidth',1.5);
    end
    ax = gca; ax.ColorOrderIndex = 1;
    for ii = 1:length(C)
        plt(2,ii) = plot(f,Pelec_f(:,ii),'-','LineWidth',1.5);
    end
    
    if jj == 1
        l1 = legend([plt(2,:)],[legCel(:)']);
    end
    set(l1,'location','northwest')
%     ylim([-5,1])
    xlim([0.2, 1])
    
    ylabel('Power [W]')
    xlabel('Frequency [Hz]')
    
    if savefigflag
        fname = sprintf('codesign_power_%s.pdf',Fe_name{jj});
        exportgraphics(gcf,fname,'ContentType','vector')
    end
    
end

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
