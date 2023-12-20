% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%%

% calculate excitation models
clear; close all;

% v3 loads wave calibration data from MASK3

%% file location information
waveBotPre ='E:\WaveBotData\waveBotCtrlMdl_2023_09_29_';
waveBotSuff = {'09_02_54','10_23_05'};
twePre = 'E:/WaveBotData/AWDC to Share/twe/awdc-t05-s000';
tweSuff = {'67','65'}; % 66 wave data set is incomplete
tdCheck =0; % to make time domain check plots

%% sys ID parameters
ssTime = [180; 30];
Trep = 300;
f_lims=[0.2 1]; % the MS excited frequency bands
N_per = 2; % number of repeat periods to grab
dt = 1/200; % slow scopes at 200 Hz
dtWave = 1/50; % wave data acquisition at 20 Hz
dtFast = 1/1000;
waveSensor = 'BUOY04';%'BRIDGEPROBE04';'BRIDGEPROBE05';...
%'BRIDGEPROBE06';'BRIDGEPROBE08';'BUOY01';

%numRuns = max(waveBotFull.expCtrl.runCounter);
Zi = importdata('./Z20A/Zikn12.mat'); % these identified systems have the same mag spring -9228
Zd = importdata('./Z20A/Zdkn12.mat'); % these identified systems have the same mag spring -9228
Zid = importdata('./Z20A/Zidkn12.mat');
load('./Z20A/WEC_T3R2.mat'); % this is BEM data for comparison

%% load wave calibration data
cd ./waveCalibration
[BPDAQ,~,~]= importNSWCCD(113); % this is a legacy function: utilize for reading MASK3 calibration ONLY!!!
cd ..

%% load data
% create data structure for each
for k=1:length(waveBotSuff)
    waveBotFull{k} = importdata([waveBotPre waveBotSuff{k} 'full.mat']);
    waveBotFull{k}.waveData.calData = BPDAQ.BRIDGEPROBE06;
    waveBotFull{k}.waveData.calBuoy = BPDAQ.(waveSensor);
    waveBotFull{k}.waveData.calWM = BPDAQ.WAVEMAKERPULSE;
    waveBotFull{k}.waveData.time = BPDAQ.Time;

    clear waveBotData;
    %% build test matrix
    for kk=1:2
        % time domain
        runSidx(kk) = find(waveBotFull{k}.expCtrl.runCounter==kk,1,'first') + ssTime(kk)/dt;
        runSidxFast(kk) = find(waveBotFull{k}.expCtrl.runCounter==kk,1,'first').*(dt/dtFast) + ssTime(kk)/dtFast;
        fpto1(k,kk,:) = waveBotFull{k}.heaveSignals.heaveForceTRS_N(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1);
        fc1(k,kk,:) = waveBotFull{k}.heaveSignals.heaveForceI_N(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1);
        ref(k,kk,:) = waveBotFull{k}.refSignals.heaveRefRaw_N(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1);
        pos(k,kk,:) = waveBotFull{k}.heaveSignals.heaveEncPos_m(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1);
        calIdx = find(waveBotFull{k}.waveData.calWM > 2.5,1,'first');
        runIdx = find(waveBotFull{k}.waveData.WAVEMAKERPULSE > -0.05,1,'first');
        if kk ==1 && k ==1; % there is only one wave calibration data set
            [~,waveT1idx] = min(abs(waveBotFull{k}.waveData.syncTimeTTL-waveBotFull{k}.heaveSignals.time(runSidx(kk))));
            [~,waveT2idx] = min(abs(waveBotFull{k}.waveData.syncTimeTTL-waveBotFull{k}.heaveSignals.time(runSidx(kk)+ N_per.*Trep/dt-1)));
                        [calR,calLags]=xcorr(waveBotFull{k}.waveData.BRIDGEPROBE06,...
                            waveBotFull{1}.waveData.calData);
                        [~,lagIdx] = max(calR);
            etaAll = interp1(waveBotFull{k}.waveData.syncTimeTTL(waveT1idx + calIdx - runIdx : waveT2idx + calIdx - runIdx),...
                waveBotFull{k}.waveData.calBuoy(waveT1idx:waveT2idx),...
                waveBotFull{k}.heaveSignals.time(runSidx(kk)+1 + (calIdx - runIdx)*(dtWave/dt) :runSidx(kk)+1 + (calIdx - runIdx)*(dtWave/dt) + N_per.*Trep/dt-1));
        end
        % time domain check plots
        if tdCheck ==1;
            figure; clf
            subplot(4,1,1)% pto
            plot(waveBotFull{k}.heaveSignals.time,waveBotFull{k}.heaveSignals.heaveForceTRS_N);
            hold on; grid on;
            plot(waveBotFull{k}.heaveSignals.time(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1),...
                waveBotFull{k}.heaveSignals.heaveForceTRS_N(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1));
            ylabel('Force N')
            subplot(4,1,2) % velocity
            plot(waveBotFull{k}.heaveSignals.time,waveBotFull{k}.heaveSignals.heaveEncPos_m);
            hold on; grid on;
            plot(waveBotFull{k}.heaveSignals.time(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1),...
                waveBotFull{k}.heaveSignals.heaveEncPos_m(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1));

            ylabel('Pos m')
)
            if kk ==1
                subplot(4,1,3) % waves
                [~,waveT1idx] = min(abs(waveBotFull{k}.waveData.syncTimeTTL-waveBotFull{k}.heaveSignals.time(runSidx(kk))));
                [~,waveT2idx] = min(abs(waveBotFull{k}.waveData.syncTimeTTL-waveBotFull{k}.heaveSignals.time(runSidx(kk)+ N_per.*Trep/dt-1)));
                plot(waveBotFull{k}.waveData.time + waveBotFull{k}.waveData.syncTimeTTL(1),waveBotFull{k}.waveData.calBuoy);
                hold on; grid on;
                plot(waveBotFull{k}.waveData.time(waveT1idx:waveT2idx)  + waveBotFull{k}.waveData.syncTimeTTL(1),waveBotFull{k}.waveData.calBuoy(waveT1idx:waveT2idx));
                plot(waveBotFull{k}.heaveSignals.time(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1),etaAll,':')
                xlabel('Time(s)')
                ylabel('eta (m)')
            end

        end

        % freqDomain
        [f_vec,~,~,F] = ampSpectra(waveBotFull{k}.heaveSignals.heaveForceTRS_N(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1),[dt; 2*dt]);
        [f_vec,~,~,Fc] = ampSpectra(waveBotFull{k}.heaveSignals.heaveForceI_N(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1),[dt; 2*dt]);
        [~,~,~,R] = ampSpectra(waveBotFull{k}.refSignals.heaveRefRaw_N(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1),[dt; 2*dt]);
        [~,~,~,P] = ampSpectra(waveBotFull{k}.heaveSignals.heaveEncPos_m(runSidx(kk):runSidx(kk) + N_per.*Trep/dt-1),[dt; 2*dt]);
        [~,~,~,ERRc] = ampSpectra(waveBotFull{k}.heaveSignals.heaveForceI_N(runSidx(kk):runSidx(kk)  + N_per.*Trep/dt-1) ...
            - waveBotFull{k}.heaveSignals.heaveForceTRS_N(runSidx(kk):runSidx(kk)  + N_per.*Trep/dt-1),[dt; 2*dt]);
        [f_fast,~,~,Ic] = ampSpectra(waveBotFull{k}.fastLogBus.heaveCurrent_A(runSidxFast(kk):runSidxFast(kk) + N_per.*Trep/dtFast-1),[dtFast; 2*dtFast]);
        [~,~,~,E] = ampSpectra(etaAll(1:end),[dt; 2*dt]);

        FPTO(k,kk,:) = F(1:N_per:end);
        FC(k,kk,:) = Fc(1:N_per:end);
        POS(k,kk,:) = P(1:N_per:end);
        ETA(k,kk,:) = E(1:N_per:end);
        REF(k,kk,:) = R(1:N_per:end);
        I(k,kk,:) = Ic(1:N_per:end);
        ERR(k,kk,:) = ERRc(1:N_per:end);
        %clear etaAll
    end
end

% find frequency ranges
f_vec2=f_vec(1:N_per:end);
f_fast2 = f_fast(1:N_per:end);
[~,f1]= min(abs(f_vec2-f_lims(1)));
[~,f2]= min(abs(f_vec2-f_lims(2)));
[~,f1fast]= min(abs(f_fast2-f_lims(1)));
[~,f2fast]= min(abs(f_fast2-f_lims(2)));
f_mult(1,1,:) = f_vec2;


%% perform calculation
% reshape to [eta1,MS1 eta1,MS2 eta2,MS1 eta2,MS2]
FPTO = reshape(FPTO,[1,4,length(FPTO(1,1,:))]);
FC = reshape(FC,[1,4,length(FPTO(1,1,:))]);
POS = reshape(POS,[1,4,length(POS(1,1,:))]);
ETA = reshape(ETA,[1,4,length(ETA(1,1,:))]);
REF = reshape(REF,[1,4,length(REF(1,1,:))]);
ERR = reshape(ERR,[1,4,length(ERR(1,1,:))]);
I = reshape(I,[1,4,length(I(1,1,:))]);
fpto = reshape(fpto1,[1,4,length(fpto1(1,1,:))]);
fc = reshape(fc1,[1,4,length(fc1(1,1,:))]);

% CHECK FULL RANGE
figure; clf;
subplot(3,1,1)
plot(f_vec2,squeeze(abs(ETA(1,1,:))));
hold on; grid on;
plot(f_vec2,squeeze(abs(ETA(1,2,:))));
plot(f_vec2,squeeze(abs(ETA(1,3,:))));
plot(f_vec2,squeeze(abs(ETA(1,4,:))));
ylabel('ETA')
subplot(3,1,2)
plot(f_vec2,squeeze(abs(FPTO(1,1,:))));
hold on; grid on;
plot(f_vec2,squeeze(abs(FPTO(1,2,:))));
plot(f_vec2,squeeze(abs(FPTO(1,3,:))));
plot(f_vec2,squeeze(abs(FPTO(1,4,:))));
ylabel('FPTO')
subplot(3,1,3)
plot(f_vec2,squeeze(abs(POS(1,1,:))));
hold on; grid on;
plot(f_vec2,squeeze(abs(POS(1,2,:))));
plot(f_vec2,squeeze(abs(POS(1,3,:))));
plot(f_vec2,squeeze(abs(POS(1,4,:))));
ylabel('POS')
xlabel('FREQ')


for k=1:length(ETA(1,1,:));
    f_exc(1,1,k) = ((POS(1,1:4,k).*Zi(1,1,k).*(2.*pi.*f_mult(1,1,k).*1i))...
        - FPTO(1,1:4,k))/ETA(1,1:4,k);
end

%% check plots
% here we discard phase b/c wave sensor not co-located with buoy
figure; clf;
plot(f_vec2(f1:f2),squeeze(abs(ETA(1,1,f1:f2))))
hold on
plot(f_vec2(f1:f2),squeeze(abs(ETA(1,2,f1:f2))))
plot(f_vec2(f1:f2),squeeze(abs(ETA(1,3,f1:f2))))
plot(f_vec2(f1:f2),squeeze(abs(ETA(1,4,f1:f2))))
xlabel('Freq (hz)')
%ylabel('abs(f_{exc}) N/m')
ylabel('abs(ETA)')

figure; clf
plot(f_vec2(f1:f2),abs(squeeze(f_exc(1,1,f1:f2))))
hold on; grid on;
plot(WEC.freq/2/pi, squeeze(abs(WEC.Exc(:,3,1))))
xlabel('Freq (Hz)')
ylabel('abs(f_{exc}) N/m')
legend('f_{SID}','f_{BEM}')
%export_fig('excEst.pdf','-transparent')
)
%% tfest system ID

H_frd = frd(squeeze(f_exc(1,1,f1:f2)),f_vec2(f1:f2),'FrequencyUnit','Hz');
H_tf = tfest(H_frd,1,0);
H_tf_frf = freqresp(H_tf,f_vec2(f1:f2),'Hz');
plot(f_vec2(f1:f2), squeeze(abs(H_tf_frf)))

%export_fig()

figure; plot(f_vec2(f1:f2), mag2db(squeeze(abs(H_tf_frf))))
hold on; grid on;
plot(f_vec2(f1:f2),mag2db(abs(squeeze(f_exc(1,1,f1:f2)))))

figure; 
subplot(2,1,1)
plot(f_vec(f1:f2),squeeze(real(f_exc(f1:f2))))
hold on; grid on;
plot(f_vec(f1:f2),squeeze(real(H_tf_frf)))
subplot(2,1,2)
plot(f_vec(f1:f2),squeeze(imag(f_exc(f1:f2))))
hold on; grid on;
plot(f_vec(f1:f2),squeeze(imag(H_tf_frf)))

