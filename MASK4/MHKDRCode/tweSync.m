% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 


function tweSynced = tweSync(tweMat, waveBotMat);

% inputs
%   tweMat: the data structure of twe file, output of tweRead.
%   waveBotMat: the data structure of the waveBot, output of loadMask4
% outputs: 
%   tweSynced: a subset of tweMat time series that is synchronized to
%       waveBotMat

wbHigh = find(waveBotMat.setpointSignals.trigger > 5,1,'first');
tweHigh = find(tweMat.TTLTrigger > 4,1,'first');
wbHight = waveBotMat.setpointSignals.time(wbHigh);
tweHight = tweMat.ElapsedTime(tweHigh);

%% this is a single-point synchronization based on the TTL signal 
% find separation in elapsed time: sync via TTL trigger
tDiff = wbHight - tweHight; 
%tweSynced = tweMat;
tweMat.syncTimeTTL = tweMat.ElapsedTime + tDiff; 

[~,sIdx] = min(abs(waveBotMat.setpointSignals.time - tweMat.syncTimeTTL(1)));
[~,eIdx] = min(abs(waveBotMat.setpointSignals.time - tweMat.syncTimeTTL(end)));

%% this is a test-averaged synchronization based on cross-correlation of the sine wave signals.
% fine adjustment using sine wave: sync via sine Wave 
SINESIGINTERP = interp1(tweMat.syncTimeTTL,tweMat.SINESIGNAL,waveBotMat.setpointSignals.time(sIdx:eIdx));
% max correlation over one period (200 samples)
[r,lags]=xcorr(SINESIGINTERP,waveBotMat.setpointSignals.syncSine(sIdx:eIdx));
[~,idx] = max(r);
if lags(idx) ~= 0;
    tweMat.syncTimeSine = tweMat.syncTimeTTL - lags(idx).*mean(diff(waveBotMat.setpointSignals.time(sIdx:eIdx)));
end

%% this is a single-point synchronization based on the GPS time
if isfield(waveBotMat.setpointSignals,'systemTime_ns')
    [~,gpsIdx] = min(abs(double(waveBotMat.setpointSignals.systemTime_ns)-tweMat.TAITime(tweHigh)));
    tGPS = waveBotMat.setpointSignals.time(gpsIdx);
    tweMat.syncTimeGPS = tweMat.ElapsedTime + (tGPS-tweMat.ElapsedTime(tweHigh));
end

%% create output structure
tweSynced = tweMat;

end