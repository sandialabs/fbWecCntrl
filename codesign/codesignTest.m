function tests = codesignTest
    clc
    tests = functiontests(localfunctions);
end

function setupOnce(testcase)
    % import matlab.unittest.fixtures.SuppressedWarningsFixture
    % testcase.applyFixture(SuppressedWarningsFixture('WAFO:MKJONSWAP'));
end

function teardownOnce(testcase)
    % Do nothing
end

%% Tests

function test_GenVsPtoImpedance(testcase)
    w = [1,2,3];
    Zgen = Gen_impedance(w,[0,1,1e-3,0]); % [Ir, Kt, Rw, Lw]
    Zpto = PTO_Impedance(w,[1, 0, 0, 0, 1, 1e-3, 0]); % [N, Id, Bd, Kd, Kt, Rw, Lw]
    verifyEqual(testcase,Zgen,Zpto)
end

function test_PowerLimit(testcase)
    cf = 60;
    mf = load('waveBot_heaveModel.mat');
    Zi = mf.Zi_frf(cf:end,1);
    Hex = mf.H_frf(cf:end,1)*1e1;
    f = mf.f(cf:end,1);
    w = 2*pi*f;
    dw = w(2)-w(1);
    
    Zpto = PTO_Impedance(w,[1, 0, 0, 0, sqrt(2/3), 1e-3, 0]);
    Fe = sqrt(8*real(Zi))*1; % const. power excitation
    
    Pmax = abs(Fe).^2 ./ (8*real(Zi));
    
    C = conj(Zi);
    ZL = Zi2ZL(Zpto,C);
    Zin = input_impedance(Zpto,ZL);
    Pmech = oneDof_mech_power(Zi, Zin, Fe);
    
    verifyEqual(testcase,Pmech,Pmax,'RelTol',1e-10)
end

function test_PiWorseThanCc_Mech(testcase)
    
    optimOpts = optimoptions('fminunc',...
    'MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');
    
    cf = 60;
    mf = load('waveBot_heaveModel.mat');
    Zi = mf.Zi_frf(cf:end,1);
    Hex = mf.H_frf(cf:end,1)*1e1;
    f = mf.f(cf:end,1);
    w = 2*pi*f;
    dw = w(2)-w(1);
    
    Zpto = PTO_Impedance(w,[1, 0, 0, 0, sqrt(2/3), 1e-3, 0]);
    Fe = sqrt(8*real(Zi))*1; % const. power excitation
    
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
    
    for ii = 1:length(wc)
        [Pmech_tot(ii), mPmech(:,ii)] = Pmech(wc(ii).ZL, Zpto, Zi, Fe);
    end
    
    verifyLessThan(testcase,-1*Pmech_tot(2),-1*Pmech_tot(1))
end

function test_PiWorseThanCc_Elec(testcase)
    
    optimOpts = optimoptions('fminunc',...
    'MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');
    
    cf = 60;
    mf = load('waveBot_heaveModel.mat');
    Zi = mf.Zi_frf(cf:end,1);
    Hex = mf.H_frf(cf:end,1)*1e1;
    f = mf.f(cf:end,1);
    w = 2*pi*f;
    dw = w(2)-w(1);
    
    Zpto = PTO_Impedance(w,[1, 0, 0, 0, sqrt(2/3), 1e-3, 0]);
    Fe = sqrt(8*real(Zi))*1; % const. power excitation
    
    %---------------------------------
    wc(1).leg = 'CC on elec';
    wc(1).ZL = conj( squeeze(Zpto(2,2,:)) ...
        - squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:)) ...
        ./ (squeeze(Zpto(1,1,:)) + Zi) );

    %---------------------------------
    wc(2).leg = 'PI on elec';
    wc(2).cinfo.type = 'PI';

    wc(2).cinfo.w = w;
    wc(2).cinfo.x0 = ones(1,2);
    wc(2).objfun = @(x) Pelec( Zi2ZL(Zpto,fbc(x,wc(2).cinfo)),...
        Zpto,...
        Zi,Fe );
    [wc(2).y, wc(2).fval] = fminunc(wc(2).objfun, wc(2).cinfo.x0, optimOpts);
    wc(2).ZL = Zi2ZL(Zpto,fbc(wc(2).y, wc(2).cinfo));
    
    for ii = 1:length(wc)
        [Pelec_tot(ii), mPelec(:,ii)] = Pelec(wc(ii).ZL, Zpto, Zi, Fe);
    end
    
    verifyLessThan(testcase,-1*Pelec_tot(2),-1*Pelec_tot(1))
    
end

function test_Mono_PiAsGoodAsCC(testcase)
    
    optimOpts = optimoptions('fminunc',...
        'MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');
    
    cf = 60;
    mf = load('waveBot_heaveModel.mat');
    Zi = mf.Zi_frf(cf:end,1);
    Hex = mf.H_frf(cf:end,1)*1e1;
    f = mf.f(cf:end,1);
    w = 2*pi*f;
    dw = w(2)-w(1);
    
    Zpto = PTO_Impedance(w,[1, 0, 0, 0, sqrt(2/3), 1e-3, 0]);
    Fe = zeros(size(Zi));
    Fe(150) = sqrt(8*real(Zi(150)))*1; % const. power excitation
    
    
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
    
    for ii = 1:length(wc)
        [Pmech_tot(ii), mPmech(:,ii)] = Pmech(wc(ii).ZL, Zpto, Zi, Fe);
    end
    
    verifyEqual(testcase,-1*Pmech_tot(2),-1*Pmech_tot(1),'RelTol',1e-4)
end

function test_Resoance_PAsGoodAsCC(testcase)
    
    optimOpts = optimoptions('fminunc',...
        'MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');
    
    cf = 60;
    mf = load('waveBot_heaveModel.mat');
    Zi = mf.Zi_frf(cf:end,1);
    Hex = mf.H_frf(cf:end,1)*1e1;
    f = mf.f(cf:end,1);
    w = 2*pi*f;
    dw = w(2)-w(1);
    
    Zpto = PTO_Impedance(w,[1, 0, 0, 0, sqrt(2/3), 1e-3, 0]);
    Fe = zeros(size(Zi));
    [~,idx] = min(abs(imag(Zi)));
    Fe(idx) = sqrt(8*real(Zi(150)))*1; % const. power excitation
    
    
    %---------------------------------
    wc(1).leg = 'CC on mech';
    wc(1).ZL = Zi2ZL(Zpto, conj(Zi));
    
    %---------------------------------
    wc(2).leg = 'P on mech';
    wc(2).cinfo.type = 'P';
    wc(2).cinfo.w = w;
    wc(2).cinfo.x0 = ones(1,2)*0.1;
    wc(2).objfun = @(x) Pmech( Zi2ZL(Zpto,fbc(x,wc(2).cinfo)),...
        Zpto,...
        Zi,Fe );
    [wc(2).y, wc(2).fval] = fminunc(wc(2).objfun, wc(2).cinfo.x0, optimOpts);
    wc(2).ZL = Zi2ZL(Zpto,fbc(wc(2).y, wc(2).cinfo));
    
    for ii = 1:length(wc)
        [Pmech_tot(ii), mPmech(:,ii)] = Pmech(wc(ii).ZL, Zpto, Zi, Fe);
    end
    
    verifyEqual(testcase,-1*Pmech_tot(2),-1*Pmech_tot(1),'RelTol',1e-5)
end

