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

function [powStudy,fh] = runPowStudy(f,Zi,Hex,S,motorSpecs,plotflag,opts)
% Designs P and PI controllers for a given WEC (specified by Zi and Hex)
% and wave spectrum (specified by S).
% 
% Args:
%   f           frequency vector
%   Zi          velocity:torque impedance (size(Zi) = [nDof,nDof,length(f)]
%   Hex         excitation (size(Hex) = [nDof, length(f)]
%   S           Wave spectrum in WAFO format
%   motorSpecs  cell array with motor specifications {Kt, R, N}, where Kt
%               is motor torque constant, R is motor winding resistance,
%               and N is gearing
%   plotFlag    set to 1 for plots
%   opts        (optional) structure to specify controller options
%     .symFlag  set to 1 to force controller to be symmetric (MIMO)
%     .diagFlag set to 1 to force diagonal controller (MIMO)

if nargin < 7
    opts.symFlag = 1;
    opts.diagFlag = 1;
end

w = 2*pi*f;
n = length(f);

nDof = size(Hex,1);

Kt = motorSpecs{1};
R = motorSpecs{2};
N = motorSpecs{3};

%% Wave spectrum and excitation

% amplitude spectrum
dw = S.w(2) - S.w(1);
ampSpect = transpose(sqrt(2*S.S*dw));

% complex excitation spectrum
Fe = Hex .* ampSpect;

%% Objective functions

% create PI controller problem
powStudy(1).Name = 'PI';
powStudy(1).cntrl = @(x) mimoPi(x,nDof,w,opts.symFlag,opts.diagFlag);
powStudy(1).x0 = mimoX0(nDof,opts.symFlag,opts.diagFlag);

% create P ("resistive damping") controller problem
powStudy(2).Name = 'P';
powStudy(2).cntrl = @(x) repmat(diag(x),[1,1,n]);
powStudy(2).x0 = zeros(nDof,1);

% populate other fields
for ii = 1:length(powStudy)
    powStudy(ii).fn = naturalRes(f,Zi);
    powStudy(ii).R = R;
    powStudy(ii).Kt = Kt;
    powStudy(ii).Ke = 2/3*Kt;
    powStudy(ii).N = N;
    powStudy(ii).f = f;
    powStudy(ii).Zi = Zi;
    powStudy(ii).Fe = Fe;
    powStudy(ii).S = S;
    powStudy(ii).powObj = @(x) WecPower(Zi,Fe,powStudy(ii).cntrl(x),Kt,R,N);
    powStudy(ii).powModel = @(C) WecPower(Zi,Fe,C,Kt,R,N);
end

%% Optimization solution

% optimization solver options
opts = optimoptions('fminunc');
opts.Display = 'off';
opts.MaxFunEvals = 5e3;
opts.TolX = 1e-8;
opts.TolFun = 1e-8;
opts.OptimalityTolerance = 1e-8;
% opts.PlotFcns = @optimplotfval;

for ii = 1:length(powStudy)
    
    powStudy(ii).optimProb.options = opts;
    powStudy(ii).optimProb.objective = powStudy(ii).powObj;
    powStudy(ii).optimProb.x0 = powStudy(ii).x0;
    powStudy(ii).optimProb.solver = 'fminunc';
    
    % solve problem
    [x,fval,ef,outp] = fminunc(powStudy(ii).optimProb);
    
    % store results
    powStudy(ii).x = x;
    powStudy(ii).fval = fval;
    powStudy(ii).ef = ef;
    powStudy(ii).outp = outp;
    [powStudy(ii).P, powStudy(ii).P_f, powStudy(ii).Pub, powStudy(ii).Pub_f] = ...
        WecPower(Zi,Fe,...
        powStudy(ii).cntrl(powStudy(ii).x),...
        Kt,powStudy(ii).R,N);
    powStudy(ii).efc = powStudy(ii).P/powStudy(ii).Pub;
    
end

%% plotting

if plotflag
    
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
    
    CntrlP = -1*tf([powStudy(1).x(1),0],[1,0]);
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
    
    for ii = 1:2
        [pl(ii)] = addGradientPatch(ax(ii),f,Fe(1,:),0.2);
    end
    
    l2 = legend(ax(2),[pl(2),Zp(:)'],...
        'Excitation, $F_e$',...
        'Theoretical lim., $Z_i^*$',...
        'PI controller, $Z_{PI}$',...
        'P controller, $Z_{P}$');
    %     [l1] = columnlegend(3, 'Excitation, $F_e$',...
    %         'Theoretical lim., $Z_i^*$',...
    %         'PI controller, $Z_{PI}$',...
    %         'P controller, $Z_{PI}$');
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

function [C] = mimoPi(x,nDof,w,symFlag,diagFlag)

w = w(:);

nFreq = length(w);

n = length(x);

kP = x(1:n/2);
kI = x(n/2+1:n);

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


function x0 = mimoX0(nDof,symFlag,diagFlag)
% 

if nDof == 1
    nVars = 2;
else
    if diagFlag
        if symFlag
            nVars = nDof;
        else
            nVars = 2*nDof;
        end
    else
        if symFlag
            nVars = 2*sum(sum(triu(ones(nDof,nDof),1) ...
                + diag([ones(nDof/2,1);zeros(nDof/2,1)])));
        else
            nVars = 2*nDof^2;
        end
    end
end

x0 = zeros(nVars,1);

end

function [P, P_f, P_ub, Pub_f]= WecPower(Zi, Fe, C, Kt, R, gear_ratio)

% Zi = [DOF,DOF,FREQ]
% Fe = [DOF, FREQ]

% electrical constant for three-phase PMS motor
Ke = Kt*2/3;

nFreq = size(Zi,3);

Pfh = zeros(nFreq,1);
Pub_f = zeros(nFreq,1);
Omega = zeros(length(R),nFreq);

for ii = 1:nFreq
    
    Omega(:,ii) = ( Zi(:,:,ii) - C(:,:,ii)) \ Fe(:,ii);
    
    Pfh(ii) = ( (Ke*gear_ratio + (R/(gear_ratio*Kt))*C(:,:,ii)) * Omega(:,ii) )'...
              * ( ((gear_ratio*Kt)\C(:,:,ii)*Omega(:,ii)) );
          
    Pub_f(ii) = -1 * real(1/8* (Fe(:,ii)' /real(Zi(:,:,ii)))*Fe(:,ii));
end
P_f = real(Pfh)*3/4;
P_ub = sum(Pub_f);
P = sum(P_f);

end


function [fn] = naturalRes(f,Z)

nDof = size(Z,1);
for ii = 1:nDof
    Zn = squeeze(Z(ii,ii,:));
    [~,indx] = min(abs(imag(Zn)));
    fn(ii) = f(indx);
end

end

function [pl] = addGradientPatch(axh,x,y,beta)

% xl = xlim(axh);
% yl = ylim(axh);
%
% [X,Y] = meshgrid(x,yl);
% Z = repmat(abs(y),[2,1]);
% n = 100;
% cmap = flipud(colormap(bone))*0.5*0.95 + 0.5;
% colormap(cmap)
% [~, pl] = contourf(axh,X,Y,Z,n,'LineStyle','none');
% p = [];

% handle for legend entry
% pl = patch(0,0,ones(1,3)*(1 - beta));
% set(pl,'EdgeColor','none');

x = x(:);
y = abs(y(:));

% n = length(x);
%
% yn = y ./ max(y);
% [~, idx] = max(y);
% xm = x(idx);

yl = ylim(axh);

% xdat = [x;flipud(x)];
% ydat = [yl(2)*ones(n,1); yl(1)*ones(n,1)];
%
% tmp1 = (1 - yn*beta) * ones(1,3);
% tmp2 = [tmp1; flipud(tmp1)];
% cdat = shiftdim(tmp2,-1);
%
% p = patch(axh,xdat,ydat,cdat);
% set(p,'FaceColor','interp','EdgeColor','none');
% uistack(p,'bottom')

b(1) = find(y>max(y)*0.5,1);
b(2) = find(y>max(y)*0.5,1,'last');

xdat = x([b(1), b(2), b(2), b(1)]);
ydat = [yl(2), yl(2), yl(1), yl(1)];

pl = patch(axh,xdat,ydat,ones(1,3)*0.95);
set(pl,'EdgeColor','none');

% plot(axh,x(b(1))*ones(2,1),yl,'k')
% plot(axh,x(b(2))*ones(2,1),yl,'k')

uistack(pl,'bottom')
set(axh,'Layer','top')

end