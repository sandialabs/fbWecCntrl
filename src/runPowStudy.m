function [powStudy,fh] = runPowStudy(f,Zi,Hex,S,options)
% runPowStudy   designs P and PI controllers for a given WEC (specified by
% Zi and Hex) and wave spectrum (specified by S).
% 
% Args:
%   f           frequency vector
%   Zi          velocity: torque impedance, size(Zi) = [nDof,nDof,length(f)]
%   Hex         excitation, size(Hex) = [nDof, length(f)]
%   S           Wave spectrum in WAFO format
%   motorSpecs  cell array with motor specifications {Kt, R, N}, where Kt
%               is motor torque constant, R is motor winding resistance,
%               and N is gearing
%   plotFlag    set to 1 for plots
%   opts        (optional) structure to specify controller options
%     .symFlag  set to 1 to force controller to be symmetric (MIMO)
%     .diagFlag set to 1 to force diagonal controller (MIMO)
%     .UpperBoundOnly set to 1 to only calculate upper bound of power
%               (complex conjugate control)
%
% Returns:
%   powStudy    struct with results
%   fh          figure handles
%
% See also mimoPi, mimoX0

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

arguments
    f (:,1) double {mustBeReal, mustBePositive, mustBeNonNan, mustBeFinite}
    Zi (:,:,:) double {mustBeNonNan, mustBeFinite}
    Hex (:,:) double {mustBeNonNan, mustBeFinite}
    S (1,1) struct
    options.Kt (:,:) double {mustBeReal, mustBeNonnegative, mustBeNonNan, mustBeFinite, checkSizes(f,Zi,Hex,options.Kt)} = eye(size(Hex,1))
    options.R (:,:) double {mustBeReal, mustBeNonnegative, mustBeNonNan, mustBeFinite} = 0*eye(size(Hex,1))
    options.N (:,:) double {mustBeReal, mustBeNonnegative, mustBeNonNan, mustBeFinite} = eye(size(Hex,1))
    options.plotFlag (1,1) {mustBeNumericOrLogical} = 0 
    options.symFlag (1,1) {mustBeNumericOrLogical} = 1
    options.diagFlag (1,1) {mustBeNumericOrLogical} = 1
    options.upperBoundOnly (1,1) {mustBeNumericOrLogical} = 0
end

%%
w = 2*pi*f;
n = length(f);

nDof = size(Hex,1);

Kt = options.Kt;
R = options.R;
N = options.N;

%% Wave spectrum and excitation

% amplitude spectrum
dw = S.w(2) - S.w(1);
ampSpect = transpose(sqrt(2*S.S*dw));

% complex excitation spectrum
Fe = Hex .* ampSpect;

%% Objective functions

% create PI controller problem
powStudy(1).Name = 'PI';
powStudy(1).nGains = 2;

% create P ("resistive damping") controller problem
powStudy(2).Name = 'P';
powStudy(2).nGains = 1;

% auto populate other fields
for ii = 1:length(powStudy)
    powStudy(ii).cntrl = @(x) mimoPi(x,...
        powStudy(ii).nGains,nDof,w,options.symFlag,options.diagFlag);
    powStudy(ii).x0 = mimoX0(powStudy(ii).nGains,...
        nDof,options.symFlag,options.diagFlag);
    powStudy(ii).fn = naturalRes(f,Zi);
    powStudy(ii).R = R;
    powStudy(ii).Kt = Kt;
    powStudy(ii).Ke = 2/3*Kt;
    powStudy(ii).N = N;
    powStudy(ii).f = f;
    powStudy(ii).Zi = Zi;
    powStudy(ii).Hex = Hex;
    powStudy(ii).Fe = Fe;
    powStudy(ii).S = S;
    powStudy(ii).ampSpect = ampSpect;
    powStudy(ii).powObj = @(x) WecPower(Zi,Fe,powStudy(ii).cntrl(x),Kt,R,N);
    powStudy(ii).powModel = @(C) WecPower(Zi,Fe,C,Kt,R,N);
    powStudy(ii).symFlag = options.symFlag;
    powStudy(ii).diagFlag = options.diagFlag;
    powStudy(ii).nDof = nDof;
end

%% Optimization solution

% optimization solver options
optim_options = optimoptions('fminunc');
optim_options.Display = 'off';
% opts.MaxFunEvals = 5e3;
% opts.TolX = 1e-8;
% opts.TolFun = 1e-8;
% opts.OptimalityTolerance = 1e-8;
% opts.PlotFcns = @optimplotfval;

for ii = 1:length(powStudy)
    
    powStudy(ii).optimProb.options = optim_options;
    powStudy(ii).optimProb.objective = powStudy(ii).powObj;
    powStudy(ii).optimProb.x0 = powStudy(ii).x0;
    powStudy(ii).optimProb.solver = 'fminunc';
    
    
    if ~options.upperBoundOnly
        
        % solve problem
        [x,fval,ef,outp] = fminunc(powStudy(ii).optimProb);
        
        % store results
        powStudy(ii).x = x;
        powStudy(ii).gainMatrix = mimoPi(powStudy(ii).x,...
            powStudy(ii).nGains,powStudy(ii).nDof,1,...
            powStudy(ii).symFlag,powStudy(ii).diagFlag);

        [powStudy(ii).P, powStudy(ii).P_f, powStudy(ii).Pub, powStudy(ii).Pub_f] = ...
            WecPower(Zi,Fe,...
            powStudy(ii).cntrl(powStudy(ii).x),...
            Kt,powStudy(ii).R,N);
        powStudy(ii).efc = powStudy(ii).P/powStudy(ii).Pub;
    else
        
        % get upper bound only
        powStudy(ii).x = NaN(size(powStudy(ii).x0));
        powStudy(ii).gainMatrix = mimoPi(powStudy(ii).x0,...
            powStudy(ii).nGains,powStudy(ii).nDof,1,...
            powStudy(ii).symFlag,powStudy(ii).diagFlag);
        fval = NaN;
        ef = NaN;
        outp = NaN;
        
        [~, ~, powStudy(ii).Pub, powStudy(ii).Pub_f] = ...
            WecPower(Zi,Fe,...
            powStudy(ii).cntrl(powStudy(ii).x0),...
            Kt,powStudy(ii).R,N);
        
        powStudy(ii).P = powStudy(ii).Pub;
        powStudy(ii).P_f = powStudy(ii).Pub_f;
    end
    
    powStudy(ii).fval = fval;
    powStudy(ii).ef = ef;
    powStudy(ii).outp = outp;
    
    powStudy(ii).efc = powStudy(ii).P/powStudy(ii).Pub;

end

%% plotting

if options.plotFlag
    
    %----------------------------------------------------------------------
    fh(1) = figure('name','Excitation',...
        'color','white',...
        'position',75+[0 0 1250 500]*0.5);
    Hex_frd = frd(shiftdim(Hex,-1),w);
    h = bodeplot(Hex_frd);
    bpopts = getoptions(h);
    bpopts.FreqUnits = 'Hz';
    bpopts.PhaseMatching = 'on';
    bpopts.PhaseMatchingFreq = 1;
    %     bpopts.PhaseWrapping = 'on';
    bpopts.Grid = 'on';
    setoptions(h,bpopts);
    setoptions(h,'XLim',[f(1), f(end)]);
    title('')
    
    allAxes = findall(fh(1),'type','axes');
    linkaxes(allAxes,'x')
    
    
    %----------------------------------------------------------------------
    fh(2) = figure('name','Impedance',...
        'color','white',...
        'position',50+[0 0 1250 500]*0.5);
    Zi_frd = frd(Zi,w);
    h = bodeplot(Zi_frd);
    bpopts = getoptions(h);
    bpopts.FreqUnits = 'Hz';
    bpopts.PhaseMatching = 'on';
    bpopts.PhaseMatchingFreq = 1;
    bpopts.PhaseWrapping = 'on';
    bpopts.Grid = 'on';
    setoptions(h,bpopts);
    setoptions(h,'XLim',[f(1), f(end)]);
    title('')
    
    %----------------------------------------------------------------------
    fh(3) = figure('name','Normalized spectral power',...
        'color','white',...
        'position',25+[0 0 1250 500]*0.5);
    hold on
    grid on
    
    Fen = abs(powStudy(1).Fe(1,:))./max(abs(powStudy(1).Fe(1,:)));
    p(1) = area(powStudy(1).f,...
        abs(powStudy(1).Fe(1,:))./max(abs(powStudy(1).Fe(1,:))));
    uistack(p(1),'bottom')
    set(gca,'Layer','top')
    %     set(p(1),'EdgeColor','none');
    p(1).FaceColor = ones(1,3)*0.95;
    p(1).DisplayName = 'Excitation, $F_e$';
    p(1).LineWidth = 0.1;
    
    p(2)  = plot(powStudy(1).f,powStudy(1).Pub_f / min(powStudy(1).Pub_f) .* Fen','b');
    %     p(2).DisplayName = sprintf('%s, %.1f\\,W','CC',abs(sum(powStudy(1).Pub_f)));
    p(2).DisplayName = 'Theoretical lim., $Z^*_i$';
    p(2).LineWidth = 2;
    
    co = {'r','m'};
    
    for ii = 1:length(powStudy)
        x = powStudy(ii).f;
        y = (powStudy(ii).P_f / min(powStudy(ii).Pub_f)) .* Fen';
        %         le = sprintf('%s, %.1f\\,W (%.0f\\%%)',...
        %             powStudy(ii).Name,abs(powStudy(ii).P), powStudy(ii).efc*100);
        le = sprintf('%s controller (%.0f\\%%)',...
            powStudy(ii).Name,powStudy(ii).efc*100);
        p(2+ii) = plot(x,y,co{ii},'DisplayName',le,'LineWidth',2);
        
    end
    
    for jj = 1:length(powStudy(1).fn)
        pfn(jj) = plot(powStudy(1).fn(jj)*ones(2,1),ylim,...
            'k--','DisplayName','WEC natural freq., $f_n$');
    end
    
    yl = ylim;
    Te = spec2char(S,5);
    pte = plot(1/Te*ones(2,1),yl,'k-',...
        'MarkerSize',3,'DisplayName','Wave energy freq., $f_e$',...
        'MarkerEdgeColor','k');
    
    l1 = legend([p,pfn(1),pte]);
    set(l1,'interpreter','latex')
    set(l1,'location','northeast')
    xlabel('Frequency [Hz]','interpreter','latex')
    ylabel('Normalized power \& excitation [ ]','interpreter','latex')
    xlim([0.2, 1.4])
    set(gca,'YTick',0:0.25:1)
    set(gca,'xscale','log')
    set(findall(fh(3),'-property','FontSize'),'FontSize',12)
    set(l1,'FontSize',10)
    ylim([0,Inf])
    
    clear pfn pte
    
    %----------------------------------------------------------------------
    fh(4) = figure('name','Impedance matching',...
        'color','white',...
        'position',[0 0 1250 500]*0.5);
    
    Zcc = conj(squeeze(Zi(1,1,:)));
    Zm = frd(Zcc,w);
    frfZcc = squeeze(freqresp(Zm,2*pi*f));
    
    CntrlPI = -1*tf([powStudy(1).x(1),powStudy(1).x(2)],[1,0]);
    frfZpi = squeeze(freqresp(CntrlPI,2*pi*f));
    
    CntrlP = -1*tf([powStudy(2).x(1),0],[1,0]);
    frfZp = squeeze(freqresp(CntrlP,2*pi*f));
    
    
    ax(1) = subplot(2,1,1);
    set(ax(1),'XScale','log')
    grid on
    hold on
    
    ax(2) = subplot(2,1,2);
    set(ax(2),'XScale','log')
    grid on
    hold on
    
    linkaxes(ax,'x')
    
    semilogx(ax(1),f,mag2db(abs(frfZcc)),'b-',...
        'LineWidth',2);
    semilogx(ax(1),f,mag2db(abs(frfZpi)),'r-',...
        'LineWidth',2);
    semilogx(ax(1),f,mag2db(abs(frfZp)),'m-',...
        'LineWidth',2);
    ylabel(ax(1),'Mag. [dB]')
    set(ax(1),'XTickLabel',[])
    
    Zp(1) = semilogx(ax(2),f,180/pi*angle(frfZcc),'b-',...
        'LineWidth',2);
    Zp(2) = semilogx(ax(2),f,180/pi*angle(frfZpi),'r-',...
        'LineWidth',2);
    Zp(3) = semilogx(ax(2),f,180/pi*angle(frfZp),'m-',...
        'LineWidth',2);
    ylabel(ax(2),'Phase [deg]')
    xlabel(ax(2),'Frequency [Hz]')
    
    ylim(ax(2),[-100,100])
    
    for ii = 1:2
        [pl(ii)] = addAreaPatch(ax(ii),f,Fe(1,:),0.5);
    end
    
    l2 = legend(ax(2),[pl(2),Zp(:)'],...
        'Excitation, $F_e$',...
        'Theoretical lim., $Z_i^*$',...
        'PI controller, $Z_{PI}$',...
        'P controller, $Z_{P}$');
    set(l2,'location','northeast',...
        'interpreter','latex')
    
    for ii = 1:2
        yl = ylim(ax(ii));
        for jj = 1
            pfn(jj) = plot(ax(ii),powStudy(1).fn(jj)*ones(2,1),yl,...
                'k--','DisplayName','WEC nat. freq., $f_n$');
        end
        
        plot(ax(ii),1/Te*ones(2,1),yl,'k-',...
            'MarkerSize',3,'DisplayName','Wave energy freq., $f_e$',...
            'MarkerEdgeColor','k');
    end
    
    title('')
    set(findall(fh(4),'-property','FontSize'),'FontSize',12)
    set(findall(fh(4),'-property','interpreter'),'interpreter','latex')
    set(l2,'FontSize',8)
    
else
    fh = [];
    
end


end


function [fn] = naturalRes(f,Z)
% naturalRes    Finds natural resonances based on impedance
%
% Finds frequency location where imaginary part of impedance is at a
% minimum
%
% Args:
%   f       frequency vector
%   Z       impedance with size(Zi) = [nDof,nDof,nFreq]
%
% Returns
%   fn      natural frequency
% 
% -------------------------------------------------------------------------

nDof = size(Z,1);
for ii = 1:nDof
    Zn = squeeze(Z(ii,ii,:));
    [~,indx] = min(abs(imag(Zn)));
    fn(ii) = f(indx);
end

end

function [pl] = addAreaPatch(axh,x,y,thresh)
% Adds an area patch in the region where y > thresh
%
% Args.
%   axh         axis handle
%   x           x data
%   y           y data
%   thresh      threshold where if y > thresh the patch will be placed
%
% -------------------------------------------------------------------------


x = x(:);
y = abs(y(:));

yl = ylim(axh);

b(1) = find(y>max(y)*thresh,1);
b(2) = find(y>max(y)*thresh,1,'last');

xdat = x([b(1), b(2), b(2), b(1)]);
ydat = [yl(2), yl(2), yl(1), yl(1)];

pl = patch(axh,xdat,ydat,ones(1,3)*0.95);
set(pl,'EdgeColor','none');

uistack(pl,'bottom')
set(axh,'Layer','top')

end

function checkSizes(f,Zi,Hex,motorSpec)
% checkSizes    Ensure that key arguments have the correct shapes
%
% -------------------------------------------------------------------------
    
if ~isequal(size(Zi,1),size(Zi,2))
    error('Argument sizes inconsistent')
end

if ~isequal(size(Zi,1),size(Hex,1))
    error('Argument sizes inconsistent')
end

for ii = 1:3
    if ~isequal(size(Zi,1),size(motorSpec,1))
        error('Argument sizes inconsistent')
    end
end

if ~isequal(size(Zi,3),size(f,1))
    error('Argument sizes inconsistent')
end
    
end
