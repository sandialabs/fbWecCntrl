% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 


%% load data
clearvars %-except cSweep; 
close all; clc;

if ~exist('cSweep')
    cSweep = importdata('E:\WaveBotData\waveBotCtrlMdl_2023_10_04_13_48_44full.mat');
    warning('cSweep already loaded in the workspace: make sure it is the right one!')
end

T_rep = 1/0.5;
endIdx =  find(cSweep.expParamSignals.heaveKpUsed == -999,1,'first');%find(cSweep.expParamSignals.heaveKpUsed==-1500,1,'last') %
startIdx= find(cSweep.expCtrl.runCounter == 1,1,'first'); %
startT = (startIdx + 200*120)/200; %(99392-(20*200))/200; % start time in s (startIdx-4001)/200;%
endT = endIdx/200; % end time in s
dT = 20; % time step between good data s
dt = 1/200; % time step slow scope
deadSamps = 5*200;
NReps =5;

nSteps = (endT - startT)./dT; % number of steps THIS SHOULD MATCH gain matrix
figure; clf;

for k=1:nSteps
    idx=ceil(startT/dt+deadSamps+(dT/dt * (k-1)):startT/dt+deadSamps+(dT/dt * (k-1)) + (NReps*T_rep)/dt);
    KpAvg(k) = round(mean(cSweep.expParamSignals.heaveKpUsed(idx)));
    KiAvg(k) = round(mean(cSweep.expParamSignals.heaveKiUsed(idx)));
    magPerc(k) = mean(cSweep.magSpringSetpointSignals.targetSpring_N_m(idx));
    P_slow =  interp1(cSweep.hbmSignals.time,cSweep.hbmSignals.heaveDriveInP,cSweep.magSpringSignals.time);
    PAvg(k) = mean(P_slow(idx));
    PAC_slow =  interp1(cSweep.hbmSignals.time,cSweep.hbmSignals.heaveDriveOutSumP,cSweep.magSpringSignals.time);
    PACAvg(k) = mean(-1.*PAC_slow(idx));
    statorForceMin(k) = min(cSweep.magSpringSignals.forceStatorA_N(idx)+cSweep.magSpringSignals.forceStatorB_N(idx));
    FslewAvg(k)=mean(cSweep.magSpringSignals.forceSlew_N(idx));
    PMech(k) = mean((detrend(cSweep.heaveSignals.heaveForceLCB_N(idx),'constant')).*cSweep.heaveSignals.heaveEncVel_m_s(idx));
    PMechTRS(k) = mean((detrend(cSweep.heaveSignals.heaveForceTRS_N(idx),'constant')).*cSweep.heaveSignals.heaveEncVel_m_s(idx)); 
    FLCBAvg(k) = mean(cSweep.heaveSignals.heaveForceLCB_N(idx));
    FLCB2(k)= prctile(cSweep.heaveSignals.heaveForceLCB_N(idx),2);
    FLCB98(k)= prctile(cSweep.heaveSignals.heaveForceLCB_N(idx),98);
    

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
subplot(8,1,1)
bar(PACAvg)
ylabel('PCAvg')
subplot(8,1,2)
bar(PMech);
ylabel('PMech')
subplot(8,1,3)
bar(KpAvg)
ylabel('KP')
subplot(8,1,4)
bar(KiAvg)
ylabel('KI')
subplot(8,1,5)
bar(magPerc)
ylabel('magSpringRate')
subplot(8,1,6)
bar(FslewAvg)
ylabel('Slew Force avg')
subplot(8,1,7)
bar(FLCBAvg)
ylabel('Force LCB avg')
subplot(8,1,8)
bar(FLCB98)
hold on;
bar(FLCB2)
ylabel('Force LCB 2-98')


magSpringIdx1=[1:20];
magSpringIdx2=[22:41];
magSpringIdx3=[43:62];
magSpringIdx4=[63:99];

figure; clf; % first spring position
subplot(1,3,1)
F111 = scatteredInterpolant(KpAvg(magSpringIdx1).',KiAvg(magSpringIdx1).',PMechTRS(magSpringIdx1).','linear','nearest');
[KpGrd111,KiGrd111] = meshgrid(unique(KpAvg(magSpringIdx1)),unique(KiAvg(magSpringIdx1)));
V111 = F111(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx1),KiAvg(magSpringIdx1),PAvg(magSpringIdx1),[],PAvg(magSpringIdx1),'filled');
contourf(KpGrd111,KiGrd111,V111)
colormap(gca,'parula')
grid on; hold on;
scatter(KpAvg(magSpringIdx1),KiAvg(magSpringIdx1),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PMechTRS (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PDC')
subplot(1,3,3)
F121 = scatteredInterpolant(KpAvg(magSpringIdx1).',KiAvg(magSpringIdx1).',PAvg(magSpringIdx1).','linear','nearest');
V121 = F121(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx1),KiAvg(magSpringIdx1),PAvg(magSpringIdx1),[],PAvg(magSpringIdx1),'filled');
contourf(KpGrd111,KiGrd111,V121)
colormap(gca,'parula')
grid on; hold on;
scatter(KpAvg(magSpringIdx1),KiAvg(magSpringIdx1),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PDC (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PDC')
title('Spring K=-9228')
subplot(1,3,2)
F122 = scatteredInterpolant(KpAvg(magSpringIdx1).',KiAvg(magSpringIdx1).',PACAvg(magSpringIdx1).','linear','nearest');
%[KpGrd121,KiGrd111] = meshgrid(unique(KpAvg(magSpringIdx1)),unique(KiAvg(magSpringIdx1)));
V122 = F122(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx1),KiAvg(magSpringIdx1),PACAvg(magSpringIdx1),[],PACAvg(magSpringIdx1),'filled');
contourf(KpGrd111,KiGrd111,V122)
grid on; hold on;
scatter(KpAvg(magSpringIdx1),KiAvg(magSpringIdx1),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PAC (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PAC')
set(gcf, 'Position', get(0, 'Screensize'));
%export_fig('powerSurface1004134844s1.pdf','-transparent') 

figure; clf; % second spring position
subplot(1,3,1)
F111 = scatteredInterpolant(KpAvg(magSpringIdx2).',KiAvg(magSpringIdx2).',PMechTRS(magSpringIdx2).','linear','nearest');
[KpGrd111,KiGrd111] = meshgrid(unique(KpAvg(magSpringIdx2)),unique(KiAvg(magSpringIdx2)));
V111 = F111(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx2),KiAvg(magSpringIdx2),PAvg(magSpringIdx2),[],PAvg(magSpringIdx2),'filled');
contourf(KpGrd111,KiGrd111,V111)
colormap(gca,'parula')
grid on; hold on;
scatter(KpAvg(magSpringIdx2),KiAvg(magSpringIdx2),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PMechTRS (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PDC')
%title('Spring Pos 1')
subplot(1,3,3)
F121 = scatteredInterpolant(KpAvg(magSpringIdx2).',KiAvg(magSpringIdx2).',PAvg(magSpringIdx2).','linear','nearest');
V121 = F121(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx2),KiAvg(magSpringIdx2),PAvg(magSpringIdx2),[],PAvg(magSpringIdx2),'filled');
contourf(KpGrd111,KiGrd111,V121)
colormap(gca,'parula')
grid on; hold on;
scatter(KpAvg(magSpringIdx2),KiAvg(magSpringIdx2),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PDC (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PDC')
title('Spring K=-8000')
subplot(1,3,2)
F122 = scatteredInterpolant(KpAvg(magSpringIdx2).',KiAvg(magSpringIdx2).',PACAvg(magSpringIdx2).','linear','nearest');
%[KpGrd121,KiGrd111] = meshgrid(unique(KpAvg(magSpringIdx2)),unique(KiAvg(magSpringIdx2)));
V122 = F122(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx2),KiAvg(magSpringIdx2),PACAvg(magSpringIdx2),[],PACAvg(magSpringIdx2),'filled');
contourf(KpGrd111,KiGrd111,V122)
grid on; hold on;
scatter(KpAvg(magSpringIdx2),KiAvg(magSpringIdx2),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PAC (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PAC')
set(gcf, 'Position', get(0, 'Screensize'));
%export_fig('powerSurface1004134844s2.pdf','-transparent') 

figure; clf; % third spring position
subplot(1,3,1)
F111 = scatteredInterpolant(KpAvg(magSpringIdx3).',KiAvg(magSpringIdx3).',PMechTRS(magSpringIdx3).','linear','nearest');
[KpGrd111,KiGrd111] = meshgrid(unique(KpAvg(magSpringIdx3)),unique(KiAvg(magSpringIdx3)));
V111 = F111(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx3),KiAvg(magSpringIdx3),PAvg(magSpringIdx3),[],PAvg(magSpringIdx3),'filled');
contourf(KpGrd111,KiGrd111,V111)
colormap(gca,'parula')
grid on; hold on;
scatter(KpAvg(magSpringIdx3),KiAvg(magSpringIdx3),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PMechTRS (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PDC')
%title('Spring Pos 1')
subplot(1,3,3)
F121 = scatteredInterpolant(KpAvg(magSpringIdx3).',KiAvg(magSpringIdx3).',PAvg(magSpringIdx3).','linear','nearest');
V121 = F121(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx3),KiAvg(magSpringIdx3),PAvg(magSpringIdx3),[],PAvg(magSpringIdx3),'filled');
contourf(KpGrd111,KiGrd111,V121)
colormap(gca,'parula')
grid on; hold on;
scatter(KpAvg(magSpringIdx3),KiAvg(magSpringIdx3),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PDC (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PDC')
title('Spring K=-6152')
subplot(1,3,2)
F122 = scatteredInterpolant(KpAvg(magSpringIdx3).',KiAvg(magSpringIdx3).',PACAvg(magSpringIdx3).','linear','nearest');
%[KpGrd121,KiGrd111] = meshgrid(unique(KpAvg(magSpringIdx3)),unique(KiAvg(magSpringIdx3)));
V122 = F122(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx3),KiAvg(magSpringIdx3),PACAvg(magSpringIdx3),[],PACAvg(magSpringIdx3),'filled');
contourf(KpGrd111,KiGrd111,V122)
grid on; hold on;
scatter(KpAvg(magSpringIdx3),KiAvg(magSpringIdx3),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PAC (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PAC')
set(gcf, 'Position', get(0, 'Screensize'));
%export_fig('powerSurface1004134844s3.pdf','-transparent') 

figure; clf; % fourth spring position
subplot(1,3,1)
F111 = scatteredInterpolant(KpAvg(magSpringIdx4).',KiAvg(magSpringIdx4).',PMechTRS(magSpringIdx4).','linear','nearest');
[KpGrd111,KiGrd111] = meshgrid(unique(KpAvg(magSpringIdx4)),unique(KiAvg(magSpringIdx4)));
V111 = F111(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx4),KiAvg(magSpringIdx4),PAvg(magSpringIdx4),[],PAvg(magSpringIdx4),'filled');
contourf(KpGrd111,KiGrd111,V111)
colormap(gca,'parula')
grid on; hold on;
scatter(KpAvg(magSpringIdx4),KiAvg(magSpringIdx4),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PMechTRS (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PDC')
%title('Spring Pos 1')
subplot(1,3,3)
F121 = scatteredInterpolant(KpAvg(magSpringIdx4).',KiAvg(magSpringIdx4).',PAvg(magSpringIdx4).','linear','nearest');
V121 = F121(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx4),KiAvg(magSpringIdx4),PAvg(magSpringIdx4),[],PAvg(magSpringIdx4),'filled');
contourf(KpGrd111,KiGrd111,V121)
colormap(gca,'parula')
grid on; hold on;
scatter(KpAvg(magSpringIdx4),KiAvg(magSpringIdx4),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PDC (W), Avg'
set(gca,'fontsize', 14) 
%zlabel('PDC')
title('Spring K=0')
subplot(1,3,2)
F122 = scatteredInterpolant(KpAvg(magSpringIdx4).',KiAvg(magSpringIdx4).',PACAvg(magSpringIdx4).','linear','nearest');
%[KpGrd121,KiGrd111] = meshgrid(unique(KpAvg(magSpringIdx4)),unique(KiAvg(magSpringIdx4)));
V122 = F122(KpGrd111,KiGrd111);
%scatter3(KpAvg(magSpringIdx4),KiAvg(magSpringIdx4),PACAvg(magSpringIdx4),[],PACAvg(magSpringIdx4),'filled');
contourf(KpGrd111,KiGrd111,V122)
grid on; hold on;
scatter(KpAvg(magSpringIdx4),KiAvg(magSpringIdx4),'ko','filled');
xlabel('Kp')
ylabel('Ki')
set(gca,'fontsize', 14) 
c= colorbar;
c.Label.String = 'PAC (W), Avg'
%zlabel('PAC')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'fontsize', 14) 
%export_fig('powerSurface1004134844s4.pdf','-transparent') 



%% remove repeated points
% repIdx = [21,22,44,45,63,64,82,83];
% KpAvg(repIdx)=[];
% KiAvg(repIdx)=[];
% magPerc(repIdx)=[];
% PAvg(repIdx)=[];
% PACAvg(repIdx)=[];


% %% Surface plots
% 
% try
%      makeSurfacePlots_draft(3,[-1000,1153.5,-9228],1,{'Kp','Ki','magK','Power_InDC'},KpAvg,KiAvg,magPerc,PAvg);
% catch
%     makeSurfacePlots_draft(3,[-1000,1153.5,-9228],1,{'Kp','Ki','magK','Power_OutAC'},KpAvg,KiAvg,magPerc,PACAvg);
% end