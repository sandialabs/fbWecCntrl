clear
clc
options = optimoptions('fminunc','MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');

cf = 60;

mf = load('waveBot_heaveModel.mat');
Zi = mf.Zi_frf(cf:end,1);
Hex = mf.H_frf(cf:end,1)*1e1;
f = mf.f(cf:end,1);
w = 2*pi*f;
dw = w(2)-w(1);


Hs = 0.125;
Tp = 5;
gamma = 1;

S = jonswap(w, [Hs, Tp, gamma]);    % Wave energy density spectrum
A = sqrt(2*dw*S.S(:));              % wave amplitude spectrum
Fe = A .* Hex(:);
% Fe1 = zeros(size(Fe));
% Fe1(end/2) = 100;
% Fe = Fe1;

Pmax = abs(Fe).^2 ./ (8*real(Zi));

%%

clc

Zpto = PTO_Impedance(w,[1, 0, 0, 0, sqrt(2/3), 0.1, 0]); % [N, Id, Bd, Kd, Kt, Rw, Lw]


%--------------------
wc(1).leg = 'CC on mech';
wc(1).ZL = Zi2ZL(Zpto, conj(Zi));

%--------------------
wc(2).leg = 'PI on mech';
wc(2).cinfo.type = 'PI';
wc(2).cinfo.w = w;
wc(2).cinfo.x0 = ones(1,2)*0.1;
wc(2).objfun = @(x) Pmech( Zi2ZL(Zpto,fbc(x,wc(2).cinfo)),...
    Zpto,...
    Zi,Fe );
[wc(2).y, wc(2).fval] = fminunc(wc(2).objfun, wc(2).cinfo.x0, options);
wc(2).ZL = Zi2ZL(Zpto,fbc(wc(2).y, wc(2).cinfo));

%--------------------
wc(3).leg = 'CC on elec';
wc(3).ZL = conj( squeeze(Zpto(2,2,:)) ...
    - squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:)) ...
    ./ (squeeze(Zpto(1,1,:)) + Zi) );

%--------------------
wc(4).leg = 'PI on elec';
wc(4).cinfo.type = 'PI';

wc(4).cinfo.w = w;
wc(4).cinfo.x0 = ones(1,2);
wc(4).objfun = @(x) Pelec( Zi2ZL(Zpto,fbc(x,wc(4).cinfo)),...
    Zpto,...
    Zi,Fe );
[wc(4).y, wc(4).fval] = fminunc(wc(4).objfun, wc(4).cinfo.x0, options);
wc(4).ZL = Zi2ZL(Zpto,fbc(wc(4).y, wc(4).cinfo));

%--------------------
wc(5).leg = 'PI on elec';
wc(5).cinfo.type = 'PI';

wc(5).cinfo.w = w;
wc(5).cinfo.x0 = ones(1,2);
wc(5).objfun = @(x) Pelec( fbc(x,wc(5).cinfo),...
    Zpto,...
    Zi,Fe );
[wc(5).y, wc(5).fval] = fminunc(wc(5).objfun, wc(5).cinfo.x0, options);
wc(5).ZL = fbc(wc(5).y,wc(5).cinfo);

for ii = 1:length(wc)
    
    % evaluate performance
    [Pmech_tot(ii), mPmech(:,ii)] = Pmech(wc(ii).ZL, Zpto, Zi, Fe);
    %     assert(-1*Pmech_tot(ii) < sum(Pmax),sprintf('''%s'' making more mechanical power than theoretical limit',wc(1).leg))
    
    [Pelec_tot(ii), mPelec(:,ii)] = Pelec(wc(ii).ZL, Zpto, Zi, Fe);
    %     assert(-1*Pelec_tot(ii) < sum(Pmax),sprintf('''%s'' making more mechanical power than theoretical limit',wc(1).leg))
end

-1 * sum(Pmax)
Pmech_tot
Pelec_tot


figure
hold on
grid on
plot(mPmech,'*')
ax = gca;
ax.ColorOrderIndex = 1;
plot(mPelec,'o')

% figure
% hold on
% grid on
% opt = bodeoptions;
% opt.FreqUnits = 'Hz';
% opt.Grid = 'on';
% for ii = 1:length(wc)
%     bodeplot(frd(wc(ii).ZL,w),opt)
%     legCel{ii} = wc(ii).leg;
% end
% legend(legCel)

function C = fbc(x,cinfo)

switch cinfo.type
    case 'PI'
        C = x(1) - 1i*x(2)./cinfo.w(:);
    case 'P'
        C = x(1);
    case 'CC'
        C = conj(cinfo.Zi);
    otherwise
        error('Invalid value for ''cinfo.type''')
end
end

function [Pmech_tot, Pmech] = Pmech(ZL, Zpto, Zi, Fe)
Zin = squeeze(Zpto(1,1,:)) - ((squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:))) ./ (squeeze(Zpto(2,2,:)) + ZL));
Pmech = 1/2 * abs( Fe ./ (Zi + Zin) ).^2 .* real(Zin);
Pmech_tot = -1 * sum(Pmech);
end

function [Pelec_tot, Pelec] = Pelec(ZL, Zpto, Zi, Fe)

Pelec = 0.5 *real(ZL) .* abs( ( squeeze(Zpto(2,1,:)) .* Fe(:)) ./ ...
    ( (( squeeze(Zpto(1,1,:)) + Zi ) .* ( squeeze(Zpto(2,2,:)) + ZL )) ...
    - squeeze(Zpto(1,2,:) .* Zpto(2,1,:)) ) ).^2;

Pelec_tot = -1 * sum(Pelec);

end

function [ZL] = Zi2ZL(Zpto, C)
    
ZL = ((squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:)) ./ ...
    ( squeeze(Zpto(1,1,:)) - C )) - squeeze(Zpto(2,2,:)));
end