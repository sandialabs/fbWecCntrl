% identifies a 1 x 1 I/O model for heave alone based on simultaneous
% multisine actuation of both drives.

clear; %close all;
msRun = importdata('E:\WaveBotData\waveBotCtrlMdl_2023_10_02_08_32_11sys.mat');
ssTime = 10; % time in seconds into each multisine signal until "steady state" (at least ramp time) 
dt = 1/200; % slow scopes at 200 Hz
numRuns = max(msRun.expCtrl.runCounter);
T_nom = 300; % ms repeat period, nominal
f_lims=[0.2 1]; % the MS excited frequency bands
N_per = 2; % number of repeat periods to grab

%% heave and force to encoder velocity
checkCorrs = 0;
checkInput =0;
k=0
for it = [1,2]%4]
    k=k+1
    % check repeat time of input to compare to nominal
    runSidx(k) = find(msRun.expCtrl.runCounter==it,1,'first') + ssTime/dt; % run start + ssTime
    acorrIn = xcorr(msRun.heaveSignals.heaveForceI_N(runSidx:end-ssTime/dt));
    [pksIn, locsIn] = findpeaks(acorrIn);
    [pkSrtIn, idx] = sort(pksIn);
    locSrtIn = locsIn(idx);
    periodSampleIn = mean([abs(locSrtIn(end)-locSrtIn(end-1)),abs(locSrtIn(end)-locSrtIn(end-2))]);
    periodTIn(k) = periodSampleIn .* dt;   

    % output signal
    clear idx;
    acorrOut = xcorr(msRun.heaveSignals.heaveEncVel_m_s(runSidx:end-ssTime/dt));
    [pksOut,locsOut] = findpeaks(acorrOut);
    [pkSrtOut, idx] = sort(pksOut);
    locSrtOut = locsOut(idx);
    periodSampleOut = mean([abs(locSrtOut(end)-locSrtOut(end-1)),abs(locSrtOut(end)-locSrtOut(end-2))]);
    periodTOut(k) = periodSampleOut .* dt;
    
    if checkCorrs == 1;
        figure; clf;
        subplot(2,1,1)
        plot(acorrIn)
        ylabel('autoCorrelation')
        hold on; grid on;
        subplot(2,1,2)
        plot(acorrOut)
        xlabel('Lags')
        ylabel('autoCorrelation')
    end

     %% set up I/O matrices tests along colum DOF along row
    % time domain
    %u(1,k,:) = msRun.fromDriveSGSignals.currentQ_A(runSidx(k):runSidx(k)) + N_per.*T_nom/dt-1);
    u1(1,k,:) = msRun.refSignals.heaveRefRaw_N(runSidx(k):runSidx(k) + N_per.*T_nom/dt-1);
    u2(1,k,:) = msRun.heaveSignals.heaveForceTRS_N(runSidx(k):runSidx(k) + N_per.*T_nom/dt-1);
    u3(1,k,:) = msRun.heaveSignals.heaveForceLCB_N(runSidx(k):runSidx(k) + N_per.*T_nom/dt-1);
    v(1,k,:) = msRun.heaveSignals.heaveEncPos_m(runSidx(k):runSidx(k) + N_per.*T_nom/dt-1);
    % frequency domain
    %[f_vec,~,~,compSpecU] = ampSpectra(msRun.heaveSignals.heaveForceI_N(runSidx(k):runSidx(k) + N_per.*T_nom/dt-1),[dt; 2*dt]);
    [f_vec,~,~,comp1SpecU] = ampSpectra(msRun.refSignals.heaveRefRaw_N(runSidx(k):runSidx(k) + N_per.*T_nom/dt-1),[dt; 2*dt]);
    [~,~,~,comp2SpecU] = ampSpectra(msRun.heaveSignals.heaveForceTRS_N(runSidx(k):runSidx(k) + N_per.*T_nom/dt-1),[dt; 2*dt]);
     [~,~,~,comp3SpecU] = ampSpectra(msRun.heaveSignals.heaveForceLCB_N(runSidx(k):runSidx(k) + N_per.*T_nom/dt-1),[dt; 2*dt]);
    
    [~,~,~,compSpecV] = ampSpectra(msRun.heaveSignals.heaveEncPos_m(runSidx(k):runSidx(k) + N_per.*T_nom/dt-1),[dt; 2*dt]);

     %U(1,k,:)= compSpecU(1:N_per:end);
     U1(1,k,:) = comp1SpecU(1:N_per:end);
     U2(1,k,:)= comp2SpecU(1:N_per:end);
     U3(1,k,:) = comp3SpecU(1:N_per:end);
     V(1,k,:)= compSpecV(1:N_per:end);
     f_vec2=f_vec(1:N_per:end);

     if checkInput ==1
         figure; clf;
         plot(f_vec,abs(comp1SpecU))
         grid on; hold on;
         xlabel('Freq (Hz)')
         ylabel('Input |U|')
     end
     
end

% find frequency ranges
[~,f1]= min(abs(f_vec2-f_lims(1))); 
[~,f2]= min(abs(f_vec2-f_lims(2)));
f_mult(1,1,:) = f_vec2;

% %% solve force I_N to enc vel
% Y=1i.*ones([1,1,length(U(1,1,:))]); % pre-allocate complex array
% for k=1:length(U(1,1,:))
%     Y(:,:,k) = U(:,:,k).'\V(:,:,k).';
% end
% Y0 = Y;
% Y= Y.*1i.*f_mult.*2*pi;

%% solve force ref to enc vel
Y1=1i.*ones([1,1,length(U1(1,1,:))]); % pre-allocate complex array
for k=1:length(U1(1,1,:))
    Y1(:,:,k) = U1(:,:,k).'\V(:,:,k).';
end
Y10 = Y1;
Y1= Y1.*1i.*f_mult.*2*pi;

%% solve force TRS to enc vel
Y2=1i.*ones([1,1,length(U2(1,1,:))]); % pre-allocate complex array
for k=1:length(U2(1,1,:))
    Y2(:,:,k) = U2(:,:,k).'\V(:,:,k).';
end
Y20 = Y2;
Y2= Y2.*1i.*f_mult.*2*pi;

%% solve force LCB to enc vel
Y3=1i.*ones([1,1,length(U3(1,1,:))]); % pre-allocate complex array
for k=1:length(U3(1,1,:))
    Y3(:,:,k) = U3(:,:,k).'\V(:,:,k).';
end
Y30 = Y3;
Y3= Y3.*1i.*f_mult.*2*pi;

%% solve for Yd from 2-port model.
Zid = 1/Y2; 
Zi = 1/Y3;
Zd = Zid./Zi;
Yd = 1/Zd;
%% plot results
% real part
figure; clf; 
subplot(2,2,1)
semilogx(f_vec2(f1:f2),squeeze(real(Y1(1,1,f1:f2))),'-b');
hold on; grid on;
%semilogx(f_vec2(f1:f2),squeeze(real(Y1(1,1,f1:f2))),'-r');
semilogx(f_vec2(f1:f2),squeeze(real(Y2(1,1,f1:f2))),'-g');
semilogx(f_vec2(f1:f2),squeeze(real(Y3(1,1,f1:f2))),'-k');
ylabel('Re(Y(1,1))')
subplot(2,2,2)
semilogx(f_vec2(f1:f2),squeeze(abs(Y1(1,1,f1:f2))),'-b');
hold on; grid on;
%semilogx(f_vec2(f1:f2),squeeze(abs(Y1(1,1,f1:f2))),'-r');
semilogx(f_vec2(f1:f2),squeeze(abs(Y2(1,1,f1:f2))),'-g');
semilogx(f_vec2(f1:f2),squeeze(abs(Y3(1,1,f1:f2))),'-k');
ylabel('abs(Y(1,1))')
subplot(2,2,3)
semilogx(f_vec2(f1:f2),squeeze(imag(Y1(1,1,f1:f2))),'-b');
hold on; grid on;
%semilogx(f_vec2(f1:f2),squeeze(imag(Y1(1,1,f1:f2))),'-r');
semilogx(f_vec2(f1:f2),squeeze(imag(Y2(1,1,f1:f2))),'-g');
semilogx(f_vec2(f1:f2),squeeze(imag(Y3(1,1,f1:f2))),'-k');
ylabel('Im(Y(1,1))')
xlabel('Freq (Hz)')
subplot(2,2,4)
semilogx(f_vec2(f1:f2),squeeze(angle(Y1(1,1,f1:f2))),'-b');
hold on; grid on;
%semilogx(f_vec2(f1:f2),squeeze(angle(Y1(1,1,f1:f2))),'-r');
semilogx(f_vec2(f1:f2),squeeze(angle(Y2(1,1,f1:f2))),'-g');
semilogx(f_vec2(f1:f2),squeeze(angle(Y3(1,1,f1:f2))),'-k');
ylabel('angle(Y(1,1))')
xlabel('Freq (Hz)')


% 
% figure; clf; 
% subplot(2,2,1)
% semilogx(f_vec2(f1:f2),squeeze(real(Y(1,1,f1:f2))));
% hold on; grid on;
% semilogx(f_vec2(f1:f2),squeeze(real(Y1(1,1,f1:f2))));
% semilogx(f_vec2(f1:f2),squeeze(real(Y2(1,1,f1:f2))));
% ylabel('Re(Y(1,1))')
% subplot(2,2,2)
% semilogx(f_vec2(f1:f2),squeeze(abs(Y(1,1,f1:f2))));
% hold on; grid on;
% semilogx(f_vec2(f1:f2),squeeze(abs(Y1(1,1,f1:f2))));
% semilogx(f_vec2(f1:f2),squeeze(abs(Y2(1,1,f1:f2))));
% ylabel('abs(Y(1,1))')
% subplot(2,2,3)
% semilogx(f_vec2(f1:f2),squeeze(imag(Y(1,1,f1:f2).*2*pi.*f_mult(f1:f2))));
% hold on; grid on;
% semilogx(f_vec2(f1:f2),squeeze(imag(Y1(1,1,f1:f2).*2*pi.*f_mult(f1:f2))));
% semilogx(f_vec2(f1:f2),squeeze(imag(Y2(1,1,f1:f2).*2*pi.*f_mult(f1:f2))));
% ylabel('Im(Y(1,1))')
% xlabel('Freq (Hz)')
% subplot(2,2,4)
% semilogx(f_vec2(f1:f2),squeeze(angle(Y(1,1,f1:f2))));
% hold on; grid on;
% semilogx(f_vec2(f1:f2),squeeze(angle(Y1(1,1,f1:f2))));
% semilogx(f_vec2(f1:f2),squeeze(angle(Y2(1,1,f1:f2))));
% ylabel('angle(Y(1,1))')
% xlabel('Freq (Hz)')


%% tfest system ID

%Y_frd = frd(squeeze(Y0(1,1,(f1:f2))),f_vec2(f1:f2),'FrequencyUnit','Hz');
Y1_frd = frd(squeeze(Y10(1,1,(f1:f2))),f_vec2(f1:f2),'FrequencyUnit','Hz');
Y2_frd = frd(squeeze(Y20(1,1,(f1:f2))),f_vec2(f1:f2),'FrequencyUnit','Hz');
Y3_frd = frd(squeeze(Y30(1,1,(f1:f2))),f_vec2(f1:f2),'FrequencyUnit','Hz');

%Y0_tf = tfest(Y_frd,2,1);
Y10_tf = tfest(Y1_frd,2,1);
Y20_tf = tfest(Y2_frd,2,1);
Y30_tf = tfest(Y3_frd,2,1);

deriv = tf([1 0],1);

%Y_tf=series(Y0_tf,deriv)
Y1_tf=series(Y10_tf,deriv)
Y2_tf=series(Y20_tf,deriv)
Y3_tf=series(Y30_tf,deriv)


%Y_tf_frf = freqresp(Y_tf,f_vec2(f1:f2),'Hz');
Y1_tf_frf = freqresp(Y1_tf,f_vec2(f1:f2),'Hz');
Y2_tf_frf = freqresp(Y2_tf,f_vec2(f1:f2),'Hz');
Y3_tf_frf = freqresp(Y3_tf,f_vec2(f1:f2),'Hz');

%% plotting 
subplot(2,2,1)
%semilogx(f_vec2(f1:f2),squeeze(real(Y_tf_frf(1,1,:))),'--b','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(real(Y1_tf_frf(1,1,:))),'--r','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(real(Y2_tf_frf(1,1,:))),'--g','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(real(Y3_tf_frf(1,1,:))),'--k','LineWidth',1.4);
subplot(2,2,2)
%semilogx(f_vec2(f1:f2),squeeze(abs(Y_tf_frf(1,1,:))),'--b','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(abs(Y1_tf_frf(1,1,:))),'--r','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(abs(Y2_tf_frf(1,1,:))),'--g','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(abs(Y3_tf_frf(1,1,:))),'--k','LineWidth',1.4);
subplot(2,2,3)
%semilogx(f_vec2(f1:f2),squeeze(imag(Y_tf_frf(1,1,:))),'--b','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(imag(Y1_tf_frf(1,1,:))),'--r','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(imag(Y2_tf_frf(1,1,:))),'--g','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(imag(Y3_tf_frf(1,1,:))),'--k','LineWidth',1.4);
subplot(2,2,4)
%semilogx(f_vec2(f1:f2),squeeze(angle(Y_tf_frf(1,1,:))),'--b','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(angle(Y1_tf_frf(1,1,:))),'--r','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(angle(Y2_tf_frf(1,1,:))),'--g','LineWidth',1.4);
semilogx(f_vec2(f1:f2),squeeze(angle(Y3_tf_frf(1,1,:))),'--k','LineWidth',1.4);
legend('Y_{ref}','Y_{TRS}','Y_{LCB}','Ytf_{ref}','Ytf_{TRS}','Ytf_{LCB}',...
'Location','SouthOutside','Orientation','Horizontal')

