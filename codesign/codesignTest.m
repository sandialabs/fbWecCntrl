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
    Zgen = Gen_impedance([0,1,1e-3,0], w); % [Ir, Kt, Rw, Lw]
    Zpto = PTO_Impedance([1, 0, 0, 0, 1, 1e-3, 0],w); % [N, Id, Bd, Kd, Kt, Rw, Lw]
    verifyEqual(testcase,Zgen,Zpto)
end

function test_PowerLimit(testcase)
    % CC matches theoretical limit of mechanical power
    
    cf = 60;
    mf = load('waveBot_heaveModel.mat');
    Zi = mf.Zi_frf(cf:end,1);
    Hex = mf.H_frf(cf:end,1)*1e1;
    f = mf.f(cf:end,1);
    w = 2*pi*f;
    dw = w(2)-w(1);
    
    Zpto = PTO_Impedance([1, 0, 0, 0, sqrt(2/3), 1e-3, 0],w);
    Fe = sqrt(8*real(Zi))*1; % const. power excitation
    
    Pmax = abs(Fe).^2 ./ (8*real(Zi));
    
    C = conj(Zi);
    ZL = Zi2ZL(Zpto,C);
    Zin = input_impedance(Zpto,ZL);
    Pmech = oneDof_mech_power(Zi, Zin, Fe);
    
    verifyEqual(testcase,Pmech,Pmax,'RelTol',1e-10)
end

function test_PiWorseThanCc_Mech(testcase)
    % PI worse than theoretical limit in sea state, mechanical power
    
    optimOpts = optimoptions('fminunc',...
        'MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');
    
    cf = 60;
    mf = load('waveBot_heaveModel.mat');
    Zi = mf.Zi_frf(cf:end,1);
    Hex = mf.H_frf(cf:end,1)*1e1;
    f = mf.f(cf:end,1);
    w = 2*pi*f;
    dw = w(2)-w(1);
    
    Zpto = PTO_Impedance([1, 0, 0, 0, sqrt(2/3), 1e-3, 0],w);
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
    % PI worse than theoretical limit in sea state, electical power
    
    optimOpts = optimoptions('fminunc',...
        'MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');
    
    cf = 60;
    mf = load('waveBot_heaveModel.mat');
    Zi = mf.Zi_frf(cf:end,1);
    Hex = mf.H_frf(cf:end,1)*1e1;
    f = mf.f(cf:end,1);
    w = 2*pi*f;
    dw = w(2)-w(1);
    
    Zpto = PTO_Impedance([1, 0, 0, 0, sqrt(2/3), 1e-3, 0],w);
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
    % PI controller matches theoretical limit in monochromatic wave
    
    optimOpts = optimoptions('fminunc',...
        'MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');
    
    cf = 60;
    mf = load('waveBot_heaveModel.mat');
    Zi = mf.Zi_frf(cf:end,1);
    Hex = mf.H_frf(cf:end,1)*1e1;
    f = mf.f(cf:end,1);
    w = 2*pi*f;
    dw = w(2)-w(1);
    
    Zpto = PTO_Impedance([1, 0, 0, 0, sqrt(2/3), 1e-3, 0],w);
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
    % P controller matches theoretical limit at resonance
    
    optimOpts = optimoptions('fminunc',...
        'MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');
    
    cf = 60;
    mf = load('waveBot_heaveModel.mat');
    Zi = mf.Zi_frf(cf:end,1);
    Hex = mf.H_frf(cf:end,1)*1e1;
    f = mf.f(cf:end,1);
    w = 2*pi*f;
    dw = w(2)-w(1);
    
    Zpto = PTO_Impedance([1, 0, 0, 0, sqrt(2/3), 1e-3, 0],w);
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

function test_OneDOF_coDesign(testcase)
    % compare w. previous results from 'OneDOF_coDesign.mlx'
    
    prev = load('OneDOF_coDesign.mat','y');
    evalc('OneDOF_coDesign'); close all
    verifyEqual(testcase, y, prev.y,'RelTol',1e-12,...
        'OneDOF_coDesign.mlx results don''t match previous')
end

function test_codesign_powerTable(testcase)
    % compare w. previous results from 'codesign_powerTable.m'
    
    evalc('codesign_powerTable'); close all;
    prev = load('codesign_powerTable.mat');
    verifyEqual(testcase, mPmech, prev.mPmech,'RelTol',1e-12,...
        'codesign_powerTable.m: mPmech results don''t match previous')
    verifyEqual(testcase, mPelec, prev.mPelec,'RelTol',1e-12,...
        'codesign_powerTable.m: mPelec results don''t match previous')
    verifyEqual(testcase, Pmax, prev.Pmax,'RelTol',1e-12,...
        'codesign_powerTable.m: mPelec results don''t match previous')
end

function test_codesign_efficiencyFig(testcase)
    % compare w. previous results from 'codesign_efficiencyFig.m'

    evalc('codesign_efficiencyFig'); close all;
    prev = load('codesign_efficiencyFig.mat');
    verifyEqual(testcase, Pelec_f, prev.Pelec_f,'RelTol',1e-12,...
        'codesign_powerTable.m: Pelec_f results don''t match previous')
    verifyEqual(testcase, Pmech_f, prev.Pmech_f,'RelTol',1e-12,...
        'codesign_powerTable.m: Pmech_f results don''t match previous')
    verifyEqual(testcase, Pmax, prev.Pmax,'RelTol',1e-12,...
        'codesign_powerTable.m: Pmax results don''t match previous')
end
