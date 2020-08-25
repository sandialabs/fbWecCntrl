clear
close all

cf = 60;

mf = load('waveBot_heaveModel.mat');
Zi = mf.Zi_frf(cf:end,1);
Hex = mf.H_frf(cf:end,1)*1e1;
f = mf.f(cf:end,1);
w = 2*pi*f;
dw = w(2)-w(1);


Hs = 0.125;
Tp = 2;
gamma = 3.3;

S = jonswap(w, [Hs, Tp, gamma]);    % Wave energy density spectrum
A = sqrt(2*dw*S.S(:));              % wave amplitude spectrum
Fe = A .* Hex(:);
% Fe = ones(size(Fe))*max(Fe);
Fe = sqrt(8*real(Zi))*1;

N = 1;  % gear ratio
Kt = sqrt(2/3); % Generator torque constant
Rw = 0; % Generator winding resistance
Lw = 0; % Generator winding inductance
Id = 0; % Drivetrain total moment of inertia (including generator's rotor)
Bd = 0; % Drivetrain friction
Kd = 0; % Drivetrain spring coefficient (between shaft and reference/ground)

PTO_params.perfect = [N, Id, Bd, Kd, Kt, Rw, Lw];

N = 12.4666;  % gear ratio
Kt = 6.1745; % Generator torque constant
Rw = 0.5; % Generator winding resistance
Lw = 1; % Generator winding inductance
Id = 200; % Drivetrain total moment of inertia (including generator's rotor)
Bd = 0.5; % Drivetrain friction
Kd = 0; % Drivetrain spring coefficient (between shaft and reference/ground)

PTO_params.actual = [N, Id, Bd, Kd, Kt, Rw, Lw];

options = optimoptions('fmincon','MaxFunctionEvaluations',1e6,'MaxIterations',1e6);

%%

Zpto = PTO_Impedance(w,[1, 0, 0, 0, 1, 1e-3, 0]); % [N, Id, Bd, Kd, Kt, Rw, Lw]

% CC on hydro
legCel{1} = 'CC on hydro';
C{1} = conj(Zi);
ZL{1} = Load_impedance(Zpto,C{1});

% CC on full sys. (eq. 22), output matching condition only
legCel{2} = 'CC on full sys.';
C{2} = conj( squeeze(Zpto(2,2,:)) ...
    - squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:)) ...
    ./ (squeeze(Zpto(1,1,:)) + Zi) );
ZL{2} = C{2};

% PI on hydro
PTO_params.ideal = [1, 0, 0, 0, 1, 1e-3, 0];
PTO_param_mask = zeros(size(PTO_params.ideal));
% 
% % generate an initial guess and bounds
x0 = ones(1,2)*0.1;


% P0 = co_optimize_oneDOF_PTO_power(x0, PTO_param_mask, PTO_params.ideal, w, Zi, Fe);
% 
coOpt_obj_fun = @(x) co_optimize_oneDOF_PTO_power(x, PTO_param_mask, PTO_params.ideal, w, Zi, Fe);
options = optimoptions('fminunc','MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6);
[y,fval] = fminunc(coOpt_obj_fun, x0, options);

Pmax = abs(Fe).^2 ./ (8*real(Zi));


for ii = 1:length(C)
    Zin{ii} = input_impedance(Zpto,ZL{ii});
    Pmech(:,ii) = oneDof_mech_power(Zi, Zin{ii}, Fe);
    [~,Pelec(:,ii)] = oneDof_PTO_power(ZL{ii},Zpto,Zi,Fe);
end



% for ii = 1:size(Pmech,2)
%     legCel{ii} = sprintf('%s (mech: %.1f W, elec: %.1f W)',...
%         legCel{ii},sum(Pmech(:,ii)),sum(Pelec(:,ii)));
% end

fig = figure;
fig.Position = fig.Position .* [1 1 1 0.5];
% set(gca,'yscale','log')
hold on
grid on

for ii = 1:length(C)
    plt(1,ii) = plot(f,Pmech(:,ii),'--','LineWidth',1.5);
end
ax = gca; ax.ColorOrderIndex = 1;
for ii = 1:length(C)
    plt(2,ii) = plot(f,Pelec(:,ii),'-','LineWidth',1.5);
end

% plot(f,Pmax,'k:')

% pd(1) = plot(0,0,'k--');
% pd(2) = plot(0,0,'k-');
% legend([pd(1),pd(2),plt(2,:)],['mech','elec', legCel(:)'])

l1 = legend([plt(2,:)],[legCel(:)']);
set(l1,'location','southwest')
ylim([-5,1])
xlim([0.2, 1])

ylabel('Efficiency [ ]')
xlabel('Frequency [Hz]')
% exportgraphics(gcf,'codesign_freqPowerComp.pdf','ContentType','vector')

%% CC on full sys.

Z_pto = PTO_Impedance(w, PTO_params.actual);
Z_L{1} = conj( Z_pto(2,2,:) - Z_pto(1,2,:) .* Z_pto(2,1,:) ./ ...
    (Z_pto(1,1,:) + shiftdim(Zi,-2)) );

%% PI on full sys.

objfun = @(x) co_optimize_oneDOF_PTO_power(x, zeros(size(PTO_params.actual)), ...
    PTO_params.actual, w, Zi, Fe);
[y,fval] = fminunc(objfun, zeros(2,1), options);
Z_L{2} = Load_impedance(w, y, Z_pto);

%% CC on hydro

Z_L{3} = conj(Zi);

%% PI on hydro

objFun = @(x) co_optimize_oneDOF_PTO_power(x, ...
    zeros(size(PTO_params)), PTO_params.perfect, w, Zi, Fe);
y = fminunc(objFun, zeros(2,1), options);
Z_L{4} = Load_impedance(w, y);

%% all powers

for ii = 1:length(Z_L)
    [P_elec(ii),~,~] = oneDof_PTO_power(squeeze(Z_L{ii}), PTO_params.actual, w, Zi, Fe);
    [P_mech(ii),~,~] = oneDof_PTO_power(squeeze(Z_L{ii}), PTO_params.perfect, w, Zi, Fe);
end

resTable = table(P_elec',P_mech',...
    'VariableNames',{'P_elec','P_mech'},...
    'RowNames',legCel);
disp(resTable)

return

%%

figure
hold on
for ii = 1:length(Z_L)
    bode(frd(Z_L{ii},w))
end
grid on
legend(legCel)
xlim([min(w),max(w)])
title('Load impedance')
