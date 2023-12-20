% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 


%%
% creates torque constant vs stator position map for mag spring from input
% stator position (X) and angular displacement (theta)
clear; close all;
%% load data

springMap = importdata('E:\WaveBotData\waveBotCtrlMdl_2023_09_19_04_44_12sys.mat');
statorPosTol = 0.01; % limit on diff(statorPosition) for a given run
filtNum = (1/20) .* ones(20,1); % this is a 10 hz moving average filter
N = 10*(0.0254); % effective radius between force @ slew and torque on magspring rotor
N_heave= 0.33047; % effective gear ratio mag spring angular rotation --> heave stage linear motion (m)
radZero = 542.5.*(pi/180);
radSignFlip= pi/4; 
deadZone = 0.1; % rads
%% step through runs
runNum = max(springMap.expCtrl.runCounter);

for k=1:runNum
%     if k==2
%         runIdx2 = find(springMap.expCtrl.runCounter==k);
%         runIdx1 = find(springMap.expCtrl.runCounter ==k-1);
%         runIdx=[runIdx1; runIdx2];
%     else
        runIdx=find(springMap.expCtrl.runCounter==k);
%     end
   
    % points when stator is stationary
    statorStats = find(diff(springMap.magSpringSignals.nanoRatePos_rev(runIdx))<statorPosTol);
    statorPos_rev(k) = median(springMap.magSpringSignals.nanoRatePos_rev(runIdx(statorStats)));
    statorPos_in(k) = median(springMap.magSpringSignals.nanoRatePos_inch(runIdx(statorStats)));
    statorPos_percent(k) = median(springMap.magSpringSignals.nanoRatePos_percent(runIdx(statorStats)));
    statorForceMax(k) = prctile(springMap.magSpringSignals.forceStatorA_N(runIdx(statorStats))+springMap.magSpringSignals.forceStatorB_N(runIdx(statorStats)),98);

    % find rotor angular displacement
    theta(:,1) = filtfilt(filtNum,1,springMap.magSpringSignals.rotorAngleEQN437_rad(runIdx(statorStats))-(springMap.magSpringSignals.slewSpAngle_deg(runIdx(statorStats)).*pi/180));
    % find torque SIGN FLIP HERE
    forceSlew(:,1) = -1.*filtfilt(filtNum,1,springMap.magSpringSignals.forceSlew_N(runIdx(statorStats)));
    torque = forceSlew .* N;
    forceHeave= torque/N_heave;
    heave_disp = theta(:,1).*N_heave;

   %% separate by displacement
   % low displacement
    hiIdx1 = find(theta < (radZero-(radSignFlip-deadZone)));
    hiIdx2 = find(theta > (radZero+(radSignFlip-deadZone)));
    hiIdx=[hiIdx1; hiIdx2];
    theta(hiIdx)=[];
    torque(hiIdx)=[];
    forceHeave(hiIdx)=[];
    heave_disp(hiIdx)=[];
    % spring Constant
    %magSpringK(k,1) = theta\torque; % N-m/rad
    %magSpringK(k,2) = theta(hiIdx)\torque(hiIdx); % N-m/rad
   
    % linear fit  FLIPPING LOAD CELL SIGN
    magSpringRad(k,[1:2]) = polyfit(theta,torque,1); % ax+b
    magSpringLin(k,[1:2]) = polyfit(heave_disp,forceHeave,1); % ax+b
    %magSpringLin(k,[3:4]) = polyfit(theta,forceSlew,1); % ax+b
   
    clear theta forceHeave forceSlew torque hiIdx hiIdx1 hiIdx2
end

figure;
plot(statorPos_percent,magSpringRad(:,1),'bo')
hold on; grid on;
xlabel('Stator Percent')
ylabel('MagSpringK(N-m/rad)')
legend('small displacement','large displacement')

figure;
plot(statorPos_percent,magSpringLin(:,1),'ro')
hold on; grid on;
xlabel('Stator Percent')
ylabel('MagSpringLin(N/m)')
legend('small displacement','large displacement')

%% format data structure for saving

magSpringMap = struct();
magSpringMap.statorPos_percent = statorPos_percent;
magSpringMap.statorPos_rev = statorPos_rev;
magSpringMap.statorPos_in = statorPos_in;
magSpringMap.magSpringKrad = magSpringRad(:,1);
magSpringMap.magSpringKLin = magSpringLin(:,1);

%magSpringMap.magSpringLin = magSpringLin; % ax+b

%%
figure; clf;
scatter3(((springMap.magSpringSignals.rotorAngleEQN437_rad(1000:end)-(springMap.magSpringSignals.slewSpAngle_deg(1000:end).*pi/180))-9.52)*180/pi,...
    springMap.magSpringSignals.forceSlew_N(1000:end).*N/N_heave,springMap.magSpringSignals.forceStatorA_N(1000:end)+springMap.magSpringSignals.forceStatorB_N(1000:end),...
    16,springMap.magSpringSignals.forceStatorA_N(1000:end)+springMap.magSpringSignals.forceStatorB_N(1000:end))
xlabel('Rot Angle deg')
ylabel('Slew force (N)')
zlabel('Stator Force (N)')
grid on
%% natural frequency calculation

mInf = 646.17;
mStat = 860;
Khyd = 24337;

wNat = sqrt((Khyd + magSpringMap.magSpringKLin)./(mInf+mStat));

figure; clf;
yyaxis left
plot(statorPos_percent,magSpringLin(:,1),'o','Color','b');
grid on; hold on;
ylabel('Spring Rate (N/m)')
yyaxis right
plot(statorPos_percent,wNat(:,1)./(2*pi),'sq','Color','r');
grid on; hold on;
xlabel('Stator Percent')
ylabel('Natural Frequency (Hz)')

%% 
