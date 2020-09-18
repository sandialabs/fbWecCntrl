% section bla bla bla
% load WEC intrinsic impedance model
clc
clear
close all

% clf
mf = load('waveBot_heaveModel.mat');
Zi = mf.Zi_frf(60:end,1);
Hex = mf.H_frf(60:end,1);
f = mf.f(60:end,1);
w = 2*pi*f;
dw = w(2)-w(1);
% Zi_frd = frd(Zi,w);
Ri = real(Zi);  % intrinsic resistance
Xi = imag(Zi);  % intrinsic reactance

% Sea state selection
% We must define a sea state on which perform the analysis
Tp = 2.5; % Wave Period (Peak period for Jonswap)
Hs = 6; % wave Height
Gamma  = 3.3; % peakiness factor (for Jonswap)


S = jonswap(w, [Hs, Tp, Gamma]);    % Wave energy density spectrum
A = sqrt(2*dw*S.S(:));              % wave amplitude spectrum
Fe = A .* Hex(:);                   % excitation spectrum

fig = figure;
fig.Position = fig.Position .* [1 1 1.25 1];


yyaxis right
area(f, 2*pi*S.S, 'facealpha',0.5)
ylabel('Spectral density [m^2/Hz]')
hold on
grid on

yyaxis left
area(f, Fe, 'facealpha',0.5)
xlim([0, 1])

ylabel('Excitation spectrum [N]')
xlabel('Frequency [Hz]')
title('Wave excitation and energy spectrum')

% PTO parameters:

% gear ratio
N = 12;
N_opt_flag = false;
N_bnds = [ 1, 50 ];
% Generator torque constant
Kt = 6.7;
Kt_opt_flag = false;
Kt_bnds = [ 3, 4 ];
% Generator winding resistance
Rw = 0.001;
Rw_opt_flag = false;
Rw_bnds = [ 3, 4 ];
% Generator winding inductance
Lw = 0;
Lw_opt_flag = false;
Lw_bnds = [ 3, 4 ];
% Drivetrain total moment of inertia (including generator's rotor)
Id = 2;
Id_opt_flag = false;
Id_bnds = [ 5, 20 ];
% Drivetrain friction
Bd = 1;
Bd_opt_flag = false;
Bd_bnds = [ 3, 4 ];
% Drivetrain spring coefficient (between shaft and reference/ground)
Kd = 0;
Kd_opt_flag = false;
Kd_bnds = [ -1e10   , 1e10 ];

PTO_cfg.PTO_param_names = {'N','Id','Bd','Kd','Kt','Rw','Lw'};
PTO_cfg.PTO_param_mask = [N_opt_flag, Id_opt_flag, Bd_opt_flag, Kd_opt_flag, Kt_opt_flag, Rw_opt_flag, Lw_opt_flag];
PTO_cfg.PTO_param_lb = [N_bnds(1), Id_bnds(1), Bd_bnds(1), Kd_bnds(1), Kt_bnds(1), Rw_bnds(1), Lw_bnds(1)];
PTO_cfg.PTO_param_ub = [N_bnds(2), Id_bnds(2), Bd_bnds(2), Kd_bnds(2), Kt_bnds(2), Rw_bnds(2), Lw_bnds(2)];

PTO_cfg.PTO_param = [N, Id, Bd, Kd, Kt, Rw, Lw];

%% no co-optimization (controller tuning only)

out_var_no_coopt = co_optimize_PTO(Fe, Zi, PTO_cfg, w);


%% co-optimization 

% Generator winding resistance
Rw = 0.1;

N_opt_flag = true;
Id_opt_flag = true;
Kd_opt_flag = false;

PTO_cfg.PTO_param_mask = [N_opt_flag, Id_opt_flag, Bd_opt_flag, Kd_opt_flag, Kt_opt_flag, Rw_opt_flag, Lw_opt_flag];
PTO_cfg.PTO_param = [N, Id, Bd, Kd, Kt, Rw, Lw];

out_var_coopt = co_optimize_PTO(Fe, Zi, PTO_cfg, w);


% ***************************************
disp(' ')
disp('******************')

fprintf('Optimal parameters for NON co-optimized PTO\n')
disp(out_var_no_coopt.T)

fprintf('Max theoretical Mechanical power [w]:\t%.2e\n',out_var_no_coopt.Pm_ub_tot)
fprintf('Max theoretical Electrical power [w]:\t%.2e\n',out_var_no_coopt.Pe_ub_tot)
fprintf('Total Electrical power [w]:\t\t%.2e\n',-1*out_var_no_coopt.P_tot)
fprintf('Theoretical max efficiency [ ]:\t\t%.2f\n', out_var_no_coopt.max_efficiency)
fprintf('Electrical efficiency [ ]:\t\t%.2f\n', out_var_no_coopt.overall_eff)
fprintf('Rel Electrical efficiency [ ]:\t\t%.2f\n', out_var_no_coopt.relative_efficiency )

disp(' ')
disp('******************')
disp(' ')

fprintf('Optimal parameters for co-optimized PTO\n')
disp(out_var_coopt.T)

fprintf('Max theoretical Mechanical power [w]:\t%.2e\n',out_var_coopt.Pm_ub_tot)
fprintf('Max theoretical Electrical power [w]:\t%.2e\n',out_var_coopt.Pe_ub_tot)
fprintf('Total Electrical power [w]:\t\t%.2e\n',-1*out_var_coopt.P_tot)
fprintf('Theoretical max efficiency [ ]:\t\t%.2f\n', out_var_coopt.max_efficiency)
fprintf('Electrical efficiency [ ]:\t\t%.2f\n', out_var_coopt.overall_eff)
fprintf('Rel Electrical efficiency [ ]:\t\t%.2f\n', out_var_coopt.relative_efficiency )

% *********************************************

function out_var = co_optimize_PTO(Fe, Zi, PTO_cfg, w)

% WEC Co-design optimization

f = w/2/pi;
Ri = real(Zi);

% generate an initial guess
x0 = [ones(1,sum(PTO_cfg.PTO_param_mask)), zeros(1,2)];

% P0 = co_optimize_oneDOF_PTO_power(x0, PTO_param_mask, PTO_param, w, Zi, Fe);

coOpt_obj_fun = @(x) co_optimize_oneDOF_PTO_power(x, PTO_cfg.PTO_param_mask, PTO_cfg.PTO_param, w, Zi, Fe);

% options = optimoptions('fminunc','MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6);
% y = fminunc(coOpt_obj_fun, x0, options)
options = optimoptions('fmincon','MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');
y = fmincon(coOpt_obj_fun, x0, [], [], [], [], ...
    PTO_cfg.PTO_param_lb(PTO_cfg.PTO_param_mask), PTO_cfg.PTO_param_ub(PTO_cfg.PTO_param_mask), [], options);

PTO_param = PTO_cfg.PTO_param;
PTO_param(PTO_cfg.PTO_param_mask) = y(1:sum(PTO_cfg.PTO_param_mask));

Z_PTO = PTO_Impedance(PTO_param, w); % PTO Impedance matrix

% plot(w, squeeze(imag(Z_PTO(1,1,:))))

% Thevenin equivalent
v_th = Fe(:) .* squeeze(Z_PTO(2,1,:)) ./ (squeeze(Z_PTO(1,1,:)) + Zi);
Z_th = squeeze(Z_PTO(2,2,:)) - squeeze((Z_PTO(1,2,:) .* Z_PTO(2,1,:) )) ...
    ./ (squeeze(Z_PTO(1,1,:) ) + Zi);

% optimal Load impedance

% plot(f, abs(v_th))
% h = bodeplot(frd(conj(Z_th),w));
% setoptions(h,'FreqUnits','Hz');
% title('Equivalent impedance Z_{th} (output impedance)')

Pm_ub = 1/8*abs(Fe).^2 ./ Ri;
out_var.Pm_ub_tot = sum(Pm_ub);

Pe_ub = 0.125 * abs(v_th(:)).^2 ./ real( Z_th );
out_var.Pe_ub_tot = sum(Pe_ub);

%  power output maching

% P_o = abs( Fe .* squeeze(Z_PTO(2,1,:)) ./ (squeeze(Z_PTO(1,1,:) ) + Zi) ).^2 ./ (8*real(Z_th));

[P_tot, P_o, ~, ~] = co_optimize_oneDOF_PTO_power(y, PTO_cfg.PTO_param_mask, PTO_param, w, Zi, Fe);

out_var.P_tot = P_tot;

eta_o = P_o ./ Pm_ub;
eta_cc = Pe_ub ./ Pm_ub;

out_var.overall_eff = -1*P_tot/out_var.Pm_ub_tot;
out_var.max_efficiency = out_var.Pe_ub_tot/out_var.Pm_ub_tot;
out_var.relative_efficiency = -1*P_tot/out_var.Pe_ub_tot;



out_var.T = table([PTO_cfg.PTO_param_lb(PTO_cfg.PTO_param_mask),nan(1,2)]',...
    [PTO_cfg.PTO_param_ub(PTO_cfg.PTO_param_mask),nan(1,2)]',...
    y',...
    'VariableNames',{'LowerBound','UpperBound','Optimal'},...
    'RowNames',[PTO_cfg.PTO_param_names(PTO_cfg.PTO_param_mask),'Kd_cntrl','Kp_cntrl']);

% fprintf('Optimal parameters\n')
% disp(T)
% 
% fprintf('Max theoretical Mechanical power [w]:\t%.2e\n',Pm_ub_tot)
% fprintf('Max theoretical Electrical power [w]:\t%.2e\n',Pe_ub_tot)
% fprintf('Total Electrical power [w]:\t\t%.2e\n',-1*P_tot)
% fprintf('Theoretical max efficiency [ ]:\t\t%.2f\n', max_efficiency)
% fprintf('Electrical efficiency [ ]:\t\t%.2f\n', overall_eff)
% fprintf('Rel Electrical efficiency [ ]:\t\t%.2f\n', relative_efficiency )
% disp(['Total power [w]: ', num2str(P_tot)])
% disp(['Overall efficiency [ ]: ', num2str(overall_eff)])



fig2 = figure;
fig2.Position = fig2.Position .* [1 1 1.25 1];

ax(1) = subplot(2,1,1);
hold on
grid on
set(ax(1),'XScale','log')
area(f, Pm_ub, 'facealpha',0.5)
area(f, Pe_ub, 'facealpha',0.5)
area(f, P_o, 'facealpha',0.5)

xlim([0.2,1])
legend(ax(1),'Theroetical Upper Bound (Mech)', 'Theroetical Upper Bound (Elec)', 'Actual Electrical power')
ylabel(ax(1),'Power [w]')

ax(2) = subplot(2,1,2);
set(ax(2),'XScale','log')
hold on
area(f, eta_cc, 'facealpha',0.5)
area(f, eta_o, 'facealpha',0.5)
ylim([0 ,1])
xlim([0.2,1])
grid on
legend(ax(2),'Theroetical Efficiency Upper Bound (CC)', 'Actual Efficiency (PI)')
ylabel(ax(2),'Efficiency [ ]')

linkaxes(ax,'x')
xlabel(ax(2),'Frequency [Hz]')

% semilogx( f, abs(den_P), f, real(den_P), f, imag(den_P))
% xlim([0.2 ,1])
% legend('abs', 'rea', 'imag')
% grid

% *************************************************************************
% *************************************************************************
end

function [P, P_f, Z_L, den_P] = co_optimize_oneDOF_PTO_power(x, PTO_param_mask, PTO_params, w, Zi, Fe)

n_PTO_var = sum(PTO_param_mask);
if n_PTO_var > 0
    PTO_params(PTO_param_mask) = x(1:n_PTO_var);
end

Z_PTO = PTO_Impedance(PTO_params, w);
Z_L = Load_impedance(x(n_PTO_var+1:end), PTO_params, w);

den_P = ((squeeze(Z_PTO(2,2,:)) + Z_L) .*...
    ( squeeze(Z_PTO(1,1,:)) + Zi )) - squeeze(Z_PTO(2,1,:).*Z_PTO(1,2,:) ) ;
P_f = 0.5 * real(Z_L) .* (abs((Fe .* squeeze(Z_PTO(2,1,:)) ) ).^2 ./ abs(den_P).^2);

P = -sum(P_f);

end

function ZL = Load_impedance(x, ~, w)

ZL = PI_cntrl(x, w);

end

function C = PI_cntrl(x, w)

C = x(1) - 1i*x(2)./w(:);

end




