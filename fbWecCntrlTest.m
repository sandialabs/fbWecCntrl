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
%% Main function to generate tests

function tests = fbWecCntrlTest
clc
tests = functiontests(localfunctions);
end

function setupOnce(testcase)
import matlab.unittest.fixtures.SuppressedWarningsFixture
testcase.applyFixture(SuppressedWarningsFixture('WAFO:MKJONSWAP'));
end

function teardownOnce(testcase)
% Do nothing
end

%% Test Functions
function test_waveBot_pow(testcase)


Hm0 = [0.127];
Tp = 1.6;
gamma = 1.0;
mf = matfile('WaveBot_heaveModel.mat');
f = mf.f;               % frequency vector
Zi = shiftdim(mf.Zi_frf,-2);
Hex = transpose(mf.H_frf);

% Create wave spectrum
Spect = jonswap(2*pi*f,[Hm0, Tp, gamma]);

Kt = 6.1745;    % WaveBot motor torque constant
R = 0.5;        % WaveBot motor electrical winding resistance
N = 12.4666;    % WaveBot heave gear ratio
motspecs = {Kt, R, N};

[powStudy] = runPowStudy(f,Zi,Hex,Spect,motspecs,0);

eval = -1.618640473816885;
verifyEqual(testcase,powStudy(1).P,eval,'RelTol',0.001)

end

function test_waveBot_gains(testcase)

Hm0 = [0.127];
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
    [powStudy] = runPowStudy(f,Zi,Hex,Spect,motspecs,0);
    
    verifyEqual(testcase,powStudy(1).x(1),evals(ii),'RelTol',0.001)
    
end
end

function test_waveBot_mono(testcase)

T = 2.5;
amp = 0.1;
H = amp*2;

Hm0 = [0.127];
Tp = [1.6];
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
[powStudy] = runPowStudy(f,Zi,Hex,Spect1,{Kt,0,N},0);

verifyEqual(testcase,powStudy(1).P,sum(powStudy(1).Pub_f),'RelTol',1e-10)

end


function test_waveBot_resonance(testcase)

T = 1/0.625;
amp = 0.1;
H = amp*2;

Hm0 = [0.127];
Tp = [1.6];
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
[powStudy] = runPowStudy(f,Zi,Hex,Spect1,{Kt,R,N},0);

verifyEqual(testcase,powStudy(1).P,powStudy(2).P,'RelTol',0.001)

end

function test_waveBot_surgePitch(testcase)

Zi = struct2array(load('Wavebot_Zi_v_enc_to_F_des_MIMO.mat'));
WEC = struct2array(load('WEC_T3R2.mat'));

Hm0 = 15*2.54e-2;
Tp = 3.5;
gamma = 1;

df = 0.01;
f = 0.2:0.01:1; % frequency vector
w = 2*pi*f;
dw = 2*pi*df;

Kt = 6.1745 * eye(2);
R = diag([0.5,  0.5]);
N = diag([12.4666, 3]);
motspecs = {Kt, R, N};

Zi_frf_sp = freqresp( Zi([2,3],[2,3]), 2*pi*f); % intrinsic impedance FRF
Zi_frf_sp(1,2,:) = Zi_frf_sp(2,1,:);

Hex = shiftdim(interp1(WEC.freq,squeeze(WEC.Exc(:,[1,5],1)),w),1);

S = jonswap(2*pi*f,[Hm0, Tp, gamma]);

opts.symFlag = 0;
opts.diagFlag = 1;

[powStudy] = runPowStudy(f,Zi_frf_sp,Hex,S,motspecs,0,opts);

eval = [-1.995366605972993, -0.204457231177428,...
    2.147752702028885, 0.130048543624792]'*1e3;

verifyEqual(testcase,powStudy(1).x,eval,'RelTol',1e-3)

end

function test_foswec_SingleFreq(testcase)

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
motspecs = {Kt, R, N};

Hm0 = 0.127;
Tp = 3;
gamma = 1.0;
S = jonswap(2*pi*f,[Hm0, Tp, gamma]);

Te = spec2char(S,5);
S2 = S;
S.S = 0*S.S;
[~,ind] = min(abs(S.w - 2*pi/Te));
S.S(ind) = S2.S(ind);

opts.symFlag = 1;
opts.diagFlag = 0;
[powStudy,~] = runPowStudy(f,Zi,Hex,S,motspecs,0,opts);

verifyEqual(testcase,powStudy(1).efc,1,'RelTol',0.01)

end

function test_foswec_gains(testcase)

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
motspecs = {Kt, R, N};

Hm0 = 0.127;
Tp = 3;
gamma = 1.0;
S = jonswap(2*pi*f,[Hm0, Tp, gamma]);

opts.symFlag = 1;
opts.diagFlag = 0;
[powStudy,~] = runPowStudy(f,Zi,Hex,S,motspecs,0,opts);

evals = [-2.276541202491610;...
    -0.923866807794542;...
    0.106576855442309;...
    -6.432190818492808];
verifyEqual(testcase,powStudy(1).x,evals,'RelTol',1e-4)

end
