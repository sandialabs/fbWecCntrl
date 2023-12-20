%% load data
clearvars %-except cSweep; 
%close all; clc;

if ~exist('cSweep')
    cSweep = loadMASK4('E:\wavebotData\waveBotCtrlMdl_2023_09_26_12_08_43.mat');
    warning('cSweep already loaded in the workspace: make sure it is the right one!')
end

T_rep = 2.5;
endIdx =  find(cSweep.expParamSignals.heaveKpUsed == -999,1,'first');%find(cSweep.expParamSignals.heaveKpUsed==-1500,1,'last') %
startIdx= find(cSweep.expCtrl.runCounter == 1,1,'first'); %
startT = (startIdx + 200*120)/200; %(99392-(20*200))/200; % start time in s (startIdx-4001)/200;%
endT = endIdx/200; % end time in s
dT = 300; % time step between good data s
dt = 1/200; % time step slow scope
deadSamps = 2*200;
NReps = 5;

nSteps = (endT - startT)./dT; % number of steps THIS SHOULD MATCH gain matrix
figure; clf;

for k=1:nSteps
    idx=ceil(startT/dt+deadSamps+(dT/dt * (k-1)):(dT/dt * (k)) + startT/dt-1);
    KpAvg(k) = round(mean(cSweep.expParamSignals.heaveKpUsed(idx)));
    KiAvg(k) = round(mean(cSweep.expParamSignals.heaveKiUsed(idx)));
    magPerc(k) = round(mean(cSweep.magSpringSetpointSignals.targetSpring_N_m(idx)));
    P_slow =  interp1(cSweep.hbmSignals.time,cSweep.hbmSignals.heaveDriveInP,cSweep.magSpringSignals.time);
    PAvg(k) = mean(P_slow(idx));
    PAC_slow =  interp1(cSweep.hbmSignals.time,cSweep.hbmSignals.heaveDriveOutSumP,cSweep.magSpringSignals.time);
    PACAvg(k) = mean(-1.*PAC_slow(idx));
    PMech(k) = mean((detrend(cSweep.heaveSignals.heaveForceTRS_N(idx),'constant')).*cSweep.heaveSignals.heaveEncVel_m_s(idx));

    if k==1;
        %plot(P_slow); 
        plot(cSweep.heaveSignals.heaveEncPos_m)
        hold on; grid on;
    else
        %plot(idx,P_slow(idx));
        plot(idx, cSweep.heaveSignals.heaveEncPos_m(idx))
    end
end

try
    makeSurfacePlots(3,[-1500,2000,-9228],1,{'Kp','Ki','magK','Power, DC (W)'},KpAvg,KiAvg,magPerc,PAvg);
catch
    makeSurfacePlots(3,[-1500,2000,-9228],1,{'Kp','Ki','magK','Power, Mech (W)'},KpAvg,KiAvg,magPerc,PMech);
end
