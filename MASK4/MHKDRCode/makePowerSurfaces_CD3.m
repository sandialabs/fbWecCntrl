% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 


%% load data
clearvars %-except cSweep; 
%close all; clc;

%% for custom drive;

if ~exist('cSweep')
    cSweep = importdata('E:\wavebotData\waveBotCtrlMdl_2023_10_03_10_21_00full.mat');
else
    warning('cSweep already loaded in the workspace: make sure it is the right one!')
end

T_rep = 1/0.3;
endIdx = find(cSweep.expParamSignals.heaveKpUsed == -999,1,'first');%find(cSweep.expCtrl.runCounter == 1,1,'last');%%find(cSweep.expParamSignals.heaveKpUsed==-1500,1,'last') %
startIdx= find(cSweep.expCtrl.runCounter == 1,1,'first'); %
startT = (startIdx + 200*120)/200; %(99392-(20*200))/200; % start time in s (startIdx-4001)/200;%
endT = endIdx/200; % end time in s
dT = 20; % time step between good data s
dt = 1/200; % time step slow scope
deadSamps = 5*200;
NReps =4;

nSteps = (endT - startT)./dT; % number of steps THIS SHOULD MATCH gain matrix
figure; clf;

dcBusVact_slow = interp1(cSweep.fromDriveSGSignals.time,cSweep.fromDriveSGSignals.voltageDC_V,cSweep.magSpringSignals.time);
PAC_slow =  interp1(cSweep.hbmSignals.time,cSweep.hbmSignals.heaveDriveOutSumP,cSweep.magSpringSignals.time);
P_slow =  interp1(cSweep.hbmSignals.time,cSweep.hbmSignals.expDriveInP,cSweep.magSpringSignals.time);
vel_slow = interp1(cSweep.fromDriveSGSignals.time,cSweep.fromDriveSGSignals.heaveEncVel_m_s,cSweep.magSpringSignals.time);    
switchT_slow = interp1(cSweep.toDriveSGSignals.time,cSweep.toDriveSGSignals.switchingPeriod_s,cSweep.magSpringSignals.time);
vel_thresh_slow = interp1(cSweep.toDriveSGSignals.time,cSweep.toDriveSGSignals.velThreshold_m_s,cSweep.magSpringSignals.time);

for k=1:nSteps
    idx=ceil(startT/dt+deadSamps+(dT/dt * (k-1)):startT/dt+deadSamps+(dT/dt * (k-1)) + (NReps*T_rep)/dt);
    KpAvg(k) = round(mean(cSweep.expParamSignals.heaveKpUsed(idx)));
    KiAvg(k) = round(mean(cSweep.expParamSignals.heaveKiUsed(idx)));
    magPerc(k) = mean(cSweep.magSpringSetpointSignals.targetSpring_N_m(idx));
    PAvg(k) = mean(P_slow(idx));
    PACAvg(k) = mean(-1.*PAC_slow(idx));
    statorForceMin(k) = min(cSweep.magSpringSignals.forceStatorA_N(idx)+cSweep.magSpringSignals.forceStatorB_N(idx));
    vel98(k)=prctile(vel_slow(idx),98);
    vel1(k)=prctile(vel_slow(idx),1);
    torque98(k)=prctile(cSweep.heaveSignals.heaveForceTRS_N(idx),98);
    torque1(k)=prctile(cSweep.heaveSignals.heaveForceTRS_N(idx),1);
    dcBusVset(k) = mean(cSweep.refSignals.dcBusVoltage_V(idx));
    dcBusVact(k) =mean(dcBusVact_slow(idx));
    driveSwitchT(k) = mean(switchT_slow(idx));
    driveVelThreshold(k) = mean(vel_thresh_slow(idx));
    
%     switch driveSwitchT(k)
%         case 0.0001
%             driveHystTol(k)=1;
%         case 0.00005
%             driveHystTol(k)=0.5;
%         case 0.00025
%             driveHystTol(k)=2;
%         case 0.0005
%             driveHystTol(k)=4;
%     end

    if k==1;
       % plot(P_slow); 
       plot(cSweep.heaveSignals.heaveEncPos_m)
       hold on; grid on;
    else
        plot(idx,cSweep.heaveSignals.heaveEncPos_m(idx));
    end
  
end



%% check plots 

figure; 
subplot(5,1,1)
bar(PAvg)
ylabel('PDCAvg (W)')
subplot(5,1,2)
bar(PACAvg)
ylabel('PACAvg (W)')
subplot(5,1,3)
bar(dcBusVset)
ylabel('VDC (V)')
subplot(5,1,4)
bar(driveSwitchT)
ylabel('driveSwitchT (s)')
subplot(5,1,5)
bar(driveVelThreshold);
ylim([0 max(driveVelThreshold)]);
ylabel('velThreshold (m/s)')
% subplot(6,1,5)
% bar(vel98)
% hold on; 
% bar(vel1)
% ylabel('Vel')
% subplot(6,1,6)
% bar(torque98)
% hold on;
% bar(torque1)
% ylabel('torque')

%chgPts = find(repPts ==1);
%chgPts = chgPts +1; 
% chgPts = [46,92,138];
% KpAvg(chgPts) = [];
% KiAvg(chgPts) = []; 
% dcBusVset(chgPts) = [];
% driveSwitchT(chgPts) = [];
% PAvg(chgPts) = [];
% PACAvg(chgPts) = [];

% magSpringIdx1=[1:30];
% magSpringIdx2=[32:54];
% magSpringIdx3=[56:77];
% magSpringIdx4=[79:103];
% 
% figure; clf; % first spring position
% subplot(1,2,1)
% scatter3(KpAvg(magSpringIdx1),KiAvg(magSpringIdx1),PAvg(magSpringIdx1),[],PAvg(magSpringIdx1),'filled');
% colormap(gca,'parula')
% grid on; hold on;
% xlabel('Kp')
% ylabel('Ki')
% zlabel('PDC')
% title('Spring Pos 1')
% subplot(1,2,2)
% scatter3(KpAvg(magSpringIdx1),KiAvg(magSpringIdx1),PACAvg(magSpringIdx1),[],PACAvg(magSpringIdx1),'filled');
% grid on; hold on;
% xlabel('Kp')
% ylabel('Ki')
% zlabel('PAC')
% 
% figure; clf; % second spring position
% subplot(1,2,1)
% scatter3(KpAvg(magSpringIdx2),KiAvg(magSpringIdx2),PAvg(magSpringIdx2),[],PAvg(magSpringIdx2),'filled');
% colormap(gca,'parula')
% grid on; hold on;
% xlabel('Kp')
% ylabel('Ki')
% zlabel('PDC')
% subplot(1,2,2)
% scatter3(KpAvg(magSpringIdx2),KiAvg(magSpringIdx2),PACAvg(magSpringIdx2),[],PACAvg(magSpringIdx2),'filled');
% grid on; hold on;
% xlabel('Kp')
% ylabel('Ki')
% zlabel('PAC')
% title('Spring Pos 2')
% 
% figure; clf; % third spring position
% subplot(1,2,1)
% scatter3(KpAvg(magSpringIdx3),KiAvg(magSpringIdx3),PAvg(magSpringIdx3),[],PAvg(magSpringIdx3),'filled');
% colormap(gca,'parula')
% grid on; hold on;
% xlabel('Kp')
% ylabel('Ki')
% zlabel('PDC')
% subplot(1,2,2)
% scatter3(KpAvg(magSpringIdx3),KiAvg(magSpringIdx3),PACAvg(magSpringIdx3),[],PACAvg(magSpringIdx3),'filled');
% grid on; hold on;
% xlabel('Kp')
% ylabel('Ki')
% zlabel('PAC')
% title('Spring Pos 3')
% 
% figure; clf; % fourth spring position
% subplot(1,2,1)
% scatter3(KpAvg(magSpringIdx4),KiAvg(magSpringIdx4),PAvg(magSpringIdx4),[],PAvg(magSpringIdx4),'filled');
% colormap(gca,'parula')
% grid on; hold on;
% xlabel('Kp')
% ylabel('Ki')
% zlabel('PDC')
% subplot(1,2,2)
% scatter3(KpAvg(magSpringIdx4),KiAvg(magSpringIdx4),PACAvg(magSpringIdx4),[],PACAvg(magSpringIdx4),'filled');
% grid on; hold on;
% xlabel('Kp')
% ylabel('Ki')
% zlabel('PAC')
% title('Spring Pos 4')



%% remove repeated points
% repIdx = [21,22,44,45,63,64,82,83];
% KpAvg(repIdx)=[];
% KiAvg(repIdx)=[];
% magPerc(repIdx)=[];
% PAvg(repIdx)=[];
% PACAvg(repIdx)=[];


% %% Surface plots

% %try
%     makeSurfacePlots_draft(3,[-2667,1538,0.1],0.01,{'Kp','Ki','velThresh','Power_InDC'},KpAvg,KiAvg,driveVelThreshold,PAvg);
% %catch
%     makeSurfacePlots_draft(3,[-2667,1538,0.1],0.01,{'Kp','Ki','velThresh','Power_OutAC'},KpAvg,KiAvg,driveVelThreshold,PACAvg);
% %end