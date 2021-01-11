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
%     along with fbWecCntrl.  If not, see <https://www.gnu.org/licenses/>.%
% -------------------------------------------------------------------------
%% Main function to generate tests

function tests = fbWecCntrlTest
clc
tests = functiontests(localfunctions);
end

function setupOnce(testcase)
import matlab.unittest.fixtures.SuppressedWarningsFixture
testcase.applyFixture(SuppressedWarningsFixture('WAFO:MKJONSWAP'));
end

function teardownOnce(~)
% Do nothing
end

%% mimoPi

function test_mimoPi_oneDof_P(testcase)

nGains = 1;
nDof = 1;
w = [1; pi];
symFlag = 0;
diagFlag = 1;
x = 12345;
C = mimoPi(x,nGains,nDof,w,symFlag,diagFlag);
wtmp = shiftdim(w,-2);
expVal = 12345 - 1i * 0 ./ wtmp;
verifyEqual(testcase,C,expVal)
    
end

function test_mimoPi_oneDof_PI(testcase)

nGains = 2;
nDof = 1;
w = [1; pi];
symFlag = 0;
diagFlag = 1;
x = [12345; 6789];
C = mimoPi(x,nGains,nDof,w,symFlag,diagFlag);
wtmp = shiftdim(w,-2);
expVal = 12345 - 1i * 6789 ./ wtmp;
verifyEqual(testcase,C,expVal)
    
end

function test_mimoPi_twoDof_P_diag(testcase)

nGains = 1;
nDof = 2;
w = [1; pi];
symFlag = 0;
diagFlag = 1;
x = [12345; 6789];
C = mimoPi(x,nGains,nDof,w,symFlag,diagFlag);
wtmp = shiftdim(w,-2);
expVal = x .* eye(2) - 1i * 0 ./ wtmp;
verifyEqual(testcase,C,expVal)
    
end

function test_mimoPi_twoDof_PI_diag(testcase)

nGains = 2;
nDof = 2;
w = [1; pi];
symFlag = 0;
diagFlag = 1;
x = [12; 34; 56; 78];
C = mimoPi(x,nGains,nDof,w,symFlag,diagFlag);
wtmp = shiftdim(w,-2);
expVal = (x(1:2) - 1i * x(3:4) ./ wtmp) .* eye(2);
verifyEqual(testcase,C,expVal)
    
end

function test_mimoPi_twoDof_P_diagSym(testcase)

nGains = 1;
nDof = 2;
w = [1; pi];
symFlag = 1;
diagFlag = 1;
x = 12345;
C = mimoPi(x,nGains,nDof,w,symFlag,diagFlag);
wtmp = shiftdim(w,-2);
expVal = x .* eye(2) - 1i * 0 ./ wtmp;
verifyEqual(testcase,C,expVal)
    
end

function test_mimoPi_twoDof_PI_diagSym(testcase)

nGains = 2;
nDof = 2;
w = [1; pi];
symFlag = 1;
diagFlag = 1;
x = [12; 34];
C = mimoPi(x,nGains,nDof,w,symFlag,diagFlag);
wtmp = shiftdim(w,-2);
expVal = (x(1) - 1i * x(2) ./ wtmp) .* eye(2);
verifyEqual(testcase,C,expVal)
    
end

function test_mimoPi_twoDof_P_sym(testcase)

nGains = 1;
nDof = 2;
w = [1; pi];
symFlag = 1;
diagFlag = 0;
x = [12345; 6789];
C = mimoPi(x,nGains,nDof,w,symFlag,diagFlag);
wtmp = shiftdim(w,-2);
expVal = [x, flipud(x)] + 1i * 0 ./ wtmp;
verifyEqual(testcase,C,expVal)
    
end

function test_mimoPi_twoDof_PI_sym(testcase)

nGains = 2;
nDof = 2;
w = [1; pi];
symFlag = 1;
diagFlag = 0;
x = [12; 34; 56; 78];
C = mimoPi(x,nGains,nDof,w,symFlag,diagFlag);
wtmp = shiftdim(w,-2);
tmp = x(1:2) - 1i * x(3:4) ./ wtmp;
expVal = [tmp, flipud(tmp)];
verifyEqual(testcase,C,expVal)
    
end

function test_mimoPi_twoDof_P_full(testcase)

nGains = 1;
nDof = 2;
w = [1; pi];
symFlag = 0;
diagFlag = 0;
x = [12; 34; 56; 78];
C = mimoPi(x,nGains,nDof,w,symFlag,diagFlag);
wtmp = shiftdim(w,-2);
expVal = reshape(x,[2,2]) - 1i * 0 ./ wtmp;
verifyEqual(testcase,C,expVal)
    
end

function test_mimoPi_twoDof_PI_full(testcase)

nGains = 2;
nDof = 2;
w = [1; pi];
symFlag = 0;
diagFlag = 0;
x = [1; 2; 3; 4; 5; 6; 7; 8];
C = mimoPi(x,nGains,nDof,w,symFlag,diagFlag);
wtmp = shiftdim(w,-2);
expVal = reshape(x(1:4),[2,2]) - 1i * reshape(x(5:8),[2,2]) ./ wtmp;
verifyEqual(testcase,C,expVal)
    
end

%% WaveBot

function test_waveBot_pow(testcase)
% test_waveBot_pow  Check solution against known power
%
% Known power established from early version of fbWecCntrl

Hm0 = 0.127;
Tp = 1.6;
gamma = 1.0;
mf = matfile('WaveBot_heaveModel.mat');
f = mf.f;               % frequency vector
Zi = shiftdim(mf.Zi_frf,-2);
Hex = transpose(mf.H_frf);

% Create wave spectrum
Spect = jonswap(2*pi*f,[Hm0, Tp, gamma]);

[powStudy] = runPowStudy(f,Zi,Hex,Spect,'Kt',6.1745,'R',0.5,'N',12.4666);

eval = -1.618640473816885;
verifyEqual(testcase,powStudy(1).P,eval,'RelTol',0.001)

end

function test_waveBot_gains(testcase)
% test_waveBot_gains    Check solution against known WaveBot heave gains
%
% Known gains established from early version of fbWecCntrl and from
% previous MASK data analysis

Hm0 = 0.127;
Tp = [1.6, 2.5, 3.5];
gamma = 1.0;
mf = matfile('WaveBot_heaveModel.mat');
f = mf.f;               % frequency vector
Zi = shiftdim(mf.Zi_frf,-2);
Hex = transpose(mf.H_frf);

evals = [-1.735081653539878, ...
    -2.280691994272000, ...
    -2.828892605607642]*1e3;

Kt = 6.1745;    % WaveBot motor torque constant
R = 0.5;        % WaveBot motor electrical winding resistance
N = 12.4666;    % WaveBot heave gear ratio
motspecs = {Kt, R, N};

for ii = 1:3
    
    % Create wave spectrum
    Spect = jonswap(2*pi*f,[Hm0, Tp(ii), gamma]);
    [powStudy] = runPowStudy(f,Zi,Hex,Spect,...
        'Kt',6.1745,'R',0.5,'N',12.4666);
    
    verifyEqual(testcase,powStudy(1).x(1),evals(ii),'RelTol',0.001)
    
end
end

function test_waveBot_mono(testcase)
% test_waveBot_mono     Ensure PI can reach optimal in monochromatic wave

T = 2.5;
amp = 0.1;

Hm0 = 0.127;
Tp = 1.6;
gamma = 1.0;
mf = matfile('WaveBot_heaveModel.mat');
f = mf.f;               % frequency vector
Zi = shiftdim(mf.Zi_frf,-2);
Hex = transpose(mf.H_frf);

Kt = 6.1745;    % WaveBot motor torque constant
R = 0.5*0;      % WaveBot motor electrical winding resistance
N = 12.4666;    % WaveBot heave gear ratio

% Create wave spectrum
Spect = jonswap(2*pi*f,[Hm0, Tp, gamma]);
[~,ind] = min(abs(Spect.w - 2*pi/T));
Spect1 = Spect;
Spect1.S = 0*Spect.S;
Spect1.w = Spect.w;
dw = Spect.w(2) - Spect.w(1);
Spect1.S(ind) = amp^2/(2*dw);
[powStudy] = runPowStudy(f,Zi,Hex,Spect1,...
        'Kt',6.1745,'R',0,'N',12.4666);

verifyEqual(testcase,powStudy(1).P,sum(powStudy(1).Pub_f),'RelTol',1e-10)

end

function test_waveBot_resonance(testcase)
% test_waveBot_resonance    At resonance (f = 0.625 Hz), P should match PI
    
T = 1/0.625;
amp = 0.1;

Hm0 = 0.127;
Tp = 1.6;
gamma = 1.0;
mf = matfile('WaveBot_heaveModel.mat');
f = mf.f;               % frequency vector
Zi = shiftdim(mf.Zi_frf,-2);
Hex = transpose(mf.H_frf);

Kt = 6.1745;    % WaveBot motor torque constant
R = 0.5;        % WaveBot motor electrical winding resistance
N = 12.4666;    % WaveBot heave gear ratio

% Create wave spectrum
Spect = jonswap(2*pi*f,[Hm0, Tp, gamma]);
[~,ind] = min(abs(Spect.w - 2*pi/T));
Spect1 = Spect;
Spect1.S = 0*Spect.S;
Spect1.w = Spect.w;
dw = Spect.w(2) - Spect.w(1);
Spect1.S(ind) = amp^2/(2*dw);
[powStudy] = runPowStudy(f,Zi,Hex,Spect1,...
        'Kt',6.1745,'R',0.5,'N',12.4666);

verifyEqual(testcase,powStudy(1).P,powStudy(2).P,'RelTol',0.001)

end

function test_waveBot_surgePitch(testcase)
% test_waveBot_surgePitch   Check solution against known WaveBot surge and
% pitch gains
%
% Known gains established from early version of fbWecCntrl and previous
% MASK data analysis

Zi = struct2array(load('waveBot_MIMO.mat'));
WEC = struct2array(load('WEC_T3R2.mat'));

Hm0 = 15*2.54e-2;
Tp = 3.5;
gamma = 1;

f = 0.2:0.01:1; % frequency vector
w = 2*pi*f;

Kt = 6.1745 * eye(2);
R = diag([0.5,  0.5]);
N = diag([12.4666, 3]);
motspecs = {Kt, R, N};

Zi_frf_sp = freqresp( Zi([2,3],[2,3]), 2*pi*f); % intrinsic impedance FRF
Zi_frf_sp(1,2,:) = Zi_frf_sp(2,1,:);

Hex = shiftdim(interp1(WEC.freq,squeeze(WEC.Exc(:,[1,5],1)),w),1);

S = jonswap(2*pi*f,[Hm0, Tp, gamma]);

[powStudy] = runPowStudy(f,Zi_frf_sp,Hex,S,...
    'symFlag',0,'diagFlag',1,'Kt',Kt,'R',R,'N',N);

eval = [-1.995366605972993, -0.204457231177428,...
    2.147752702028885, 0.130048543624792]'*1e3;

verifyEqual(testcase,powStudy(1).x,eval,'RelTol',1e-3)

end

% function test_wavebot_defaults(testcase)
% % test_wavebot_defaults     Check that default parameters work
% 
% mf = matfile('foswec_model.mat');
% f = mf.f;
% Hex = mf.Hex;
% Zi = mf.Zi;
% 
% Hm0 = 0.127;
% Tp = 3;
% gamma = 1.0;
% S = jonswap(2*pi*f,[Hm0, Tp, gamma]);
% 
% [powStudy,~] = runPowStudy(f,Zi,Hex,S,...
%     'symFlag',1,'diagFlag',1,'Kt',eye(size(Hex,1)),'R',0*eye(size(Hex,1)),'N',eye(size(Hex,1)));
% [powStudy1,~] = runPowStudy(f,Zi,Hex,S);
% 
% verifyEqual(testcase,powStudy1(1).P,powStudy(1).P,'RelTol',1e-4)
% 
% end


%% FOSWEC

function test_foswec_SingleFreq(testcase)
% test_foswec_SingleFreq    Ensure PI can reach optimal in monochromatic
% wave

mf = matfile('foswec_model.mat');
f = mf.f;
Hex = mf.Hex;
Zi = mf.Zi;

Kt = 0.943;
R = 1.082;
N = 3.75;
R = eye(size(Hex,1))*R*0;
Kt = eye(size(Hex,1))*Kt;
N = eye(size(Hex,1))*N;

Hm0 = 0.127;
Tp = 3;
gamma = 1.0;
S = jonswap(2*pi*f,[Hm0, Tp, gamma]);

Te = spec2char(S,5);
S2 = S;
S.S = 0*S.S;
[~,ind] = min(abs(S.w - 2*pi/Te));
S.S(ind) = S2.S(ind);

[powStudy,~] = runPowStudy(f,Zi,Hex,S,...
    'symFlag',1,'diagFlag',0,'Kt',Kt,'R',R,'N',N);

verifyEqual(testcase,powStudy(1).efc,1,'RelTol',0.01)

end

function test_foswec_gains(testcase)
% test_foswec_gains     Check solution against known FOSWEC gains
%
% Known gains established from early version of fbWecCntrl

mf = matfile('foswec_model.mat');
f = mf.f;
Hex = mf.Hex;
Zi = mf.Zi;

Kt = 0.943;
R = 1.082;
N = 3.75;
R = eye(size(Hex,1))*R*0;
Kt = eye(size(Hex,1))*Kt;
N = eye(size(Hex,1))*N;

Hm0 = 0.127;
Tp = 3;
gamma = 1.0;
S = jonswap(2*pi*f,[Hm0, Tp, gamma]);

[powStudy,~] = runPowStudy(f,Zi,Hex,S,...
    'symFlag',1,'diagFlag',0,'Kt',Kt,'R',R,'N',N);

expvals = [-23.807678927193628;...
    -12.234542384574750;...
    -4.641103810241376;...
    -81.659854004619802];
verifyEqual(testcase,powStudy(1).x,expvals,'RelTol',1e-4)

end


function test_foswec_defaults(testcase)
% test_foswec_defaults     Check that default parameters work

mf = matfile('foswec_model.mat');
f = mf.f;
Hex = mf.Hex;
Zi = mf.Zi;

Hm0 = 0.127;
Tp = 3;
gamma = 1.0;
S = jonswap(2*pi*f,[Hm0, Tp, gamma]);

[powStudy,~] = runPowStudy(f,Zi,Hex,S,...
    'symFlag',1,'diagFlag',1,'Kt',eye(size(Hex,1)),'R',0*eye(size(Hex,1)),'N',eye(size(Hex,1)));
[powStudy1,~] = runPowStudy(f,Zi,Hex,S);

verifyEqual(testcase,powStudy1(1).P,powStudy(1).P,'RelTol',1e-4)

end

