clear
close all

mf = load('waveBot_heaveModel.mat');
Zi = mf.Zi_frf(60:end,1);
Hex = mf.H_frf(60:end,1)*1e1;
f = mf.f(60:end,1);
w = 2*pi*f;
dw = w(2)-w(1);

Hs = 0.125;
Tp = 2;
gamma = 3.3;

S = jonswap(w, [Hs, Tp, gamma]);    % Wave energy density spectrum
A = sqrt(2*dw*S.S(:));              % wave amplitude spectrum
Fe = A .* Hex(:);

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

legCel = {...
    'CC on full sys.',...
    'PI on full sys.',...
    'CC on hydro',...
    'PI on hydro',...
    };

%% CC on full sys.

Z_pto = PTO_Impedance(w, PTO_params.actual);
Z_L{1} = conj( Z_pto(2,2,:) - Z_pto(1,2,:) .* Z_pto(2,1,:) ./ (Z_pto(1,1,:) + shiftdim(Zi,-2)) );

%% PI on full sys.

objfun = @(x) co_optimize_oneDOF_PTO_power(x, zeros(size(PTO_params.actual)), ...
    PTO_params.actual, w, Zi, Fe);
options = optimoptions('fmincon','MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6);
[y,fval] = fminunc(objfun, zeros(2,1), options);
P(2) = fval;
Z_L{2} = Load_impedance(w, y);

%% CC on hydro

Z_L{3} = conj(Zi);

%% PI on hydro

coOpt_obj_fun = @(x) co_optimize_oneDOF_PTO_power(x, zeros(size(PTO_params)), PTO_params.perfect, w, Zi, Fe);
options = optimoptions('fminunc','MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6);
y = fminunc(coOpt_obj_fun, zeros(2,1), options);
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
