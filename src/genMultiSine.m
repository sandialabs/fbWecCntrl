function [t,u,fFreq,uFreq,phFreq,p2p] = genMultiSine(fmin,fmax,T,options)
    % genMultiSine   generate multisine signal
    %
    % Generates a white or pink multi-sine signal to be used in system
    % identification (SID) experiments. A large set of phase realizations
    % are generated and the realization with the smallest peak-to-peak
    % range is returned.
    %
    % Arguments:
    %   fmin      minimum frequency [Hz]
    %   fmax      maximum frequency [Hz]
    %   T         period of signal [s]
    %
    % optional name-value pairs
    %   nPh       number of phasings to consider (default: 1e2)
    %   plotFlag  set to 1 for plots (default: 0)
    %   color     'pink' or 'white' (default: 'white')
    %   dt        time-step for time series [s] (default: 1e-2)
    %
    % Returns:
    %   t         time series vector [s]
    %   u         signal time vector [*]
    %   fFreq     frequency vector [Hz]
    %   uFreq     signal frequency vector [*]
    %   phFreq    phase vector [rad]
    %   p2p       peak-to-peak signal range [*]
    %
    % Example:
    %
    % [t,u] = genMultiSine(0.25,1,300,'plotFlag',true,'color','white');
    %
    % Authors: Giorgio Bacelli (original), Ryan Coe (minor tweaks)
    
    % ---------------------------------------------------------------------
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
        fmin (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
        fmax (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
        T (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
        options.numPhases double {mustBeFinite,mustBeReal,mustBePositive}  = 1e2
        options.plotFlag {mustBeNumericOrLogical} = 0
        options.color string = 'white'
        options.dt (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 1e-2
        options.RMS (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 1
    end
    
    df = 1/T;
    dt = options.dt;
    
    ind_fmin = round(fmin / df);
    ind_fmax = round(fmax / df);
    
    f_vec = (ind_fmin:ind_fmax)'*df;
    N_freq = length(f_vec);
    
    ph_mat = 2*pi*rand(N_freq, options.numPhases);
    
%     Amp = ones(size(f_vec));
    
    switch options.color
        case 'pink'
            Amp = ones(size(f_vec)) ./ sqrt(f_vec);
        case 'white'
            Amp = ones(size(f_vec));
        otherwise
            error('Unkown noise color (choose ''pink'' or ''white'')')
    end

    RMS_tmp = sqrt(sum(Amp.^2))/sqrt(2)/options.RMS;
    Amp = Amp / RMS_tmp;
    
    t_vec = (0:dt:T-dt)';
    exp_mat = exp(1i * 2*pi*f_vec * (t_vec'));
    
%     fd_ms = [];
%     p2p = Inf;
    
    % TODO: use this approach if the problem is really big
%     for ind_ph = 1: options.numPhases
%         fd_ms_tmp(:,ind_ph) = Amp .* exp(1i .* ph_mat(:,ind_ph));
%         td_ms_tmp = real(fd_ms_tmp.' * exp_mat);
%         delta_amp = max(td_ms_tmp) - min(td_ms_tmp);
%         if delta_amp < p2p
%             p2p = delta_amp;
%             ph = ph_mat(:, ind_ph);
%             fd_ms = fd_ms_tmp; %/max(abs(sig_tmp));
%             td_ms = td_ms_tmp;
%         end
%     end
    
    fd_ms_tmp = Amp .* exp(1i .* ph_mat);
    td_ms_tmp = real(fd_ms_tmp.' * exp_mat)';
    delta_amp = max(td_ms_tmp,[],1) - min(td_ms_tmp,[],1);
    
    [p2p,midx] = min(delta_amp);
    td_ms = td_ms_tmp(:,midx);% / (0.5*p2p);
    fd_ms = fd_ms_tmp(:,midx);
    ph_ms = ph_mat(:,midx);
    
    if options.plotFlag
        
        figure
        
        subplot 211
        loglog(f_vec, Amp, '.')
%         ylim([0,1])
        xlabel('Frequency [Hz]','interpreter','latex');
        ylabel('Amplitude','interpreter','latex');
        title(sprintf('%s noise, $f \\in [%.2g,%.2g]$, RMS = %g',...
            options.color,fmin,fmax, options.RMS),'interpreter','latex');
        grid on
        
        subplot 212
        
        plot(t_vec, td_ms)
        xlabel('Time [s]','interpreter','latex');
        ylabel('Ampiltude','interpreter','latex');
        title(sprintf('period: %3g s',...
            T),'interpreter','latex');
        grid on
%         ylim([-1.25, 1.25])
    end
    
    t = t_vec;
    u = td_ms;
    
    fFreq = f_vec;
    uFreq = fd_ms;
    phFreq = ph_ms;
end
