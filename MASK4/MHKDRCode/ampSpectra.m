% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

function [freq,ampSpec,phaseSpec,compSpec]=ampSpectra(Data,time)

%% INPUTS:
% Data: array of column vectors containing the data of interest. Ensure it
% is periodic or you'll have spectral leakage, zero pad to reduce
% time: the time series (in s) related to the Data
%% OUTPUTS:
% freq: the vector of frequencies (Hz) for the one-sided amplitude spectrum
% ampSpec: the one-sided amplitude spectrum
% phaseSpec: the phase of the one-sided spectrum
% compSpec: the complex number representation of the 

% detrend, zero pad, and fft
%Data=detrend(Data);
L=length(Data);
%Ln=2^nextpow2(L);
Ln=L;
Y=fft(Data,Ln);

% handles frequency vector
dt=mean(diff(time));
if var(diff(time)) > 1e-5;
    warning('Non-uniform sample time. Must resample data and try again')
end
Fs=1/dt;
freq=Fs*(0:(Ln/2))/Ln;
%freq=freq(2:end-1);

% calculated amplitude spectra
P2=(Y/Ln);
P1=P2(1:(Ln/2)+1);
P1(2:end-1)=2*P1(2:end-1);
compSpec = P1;
ampSpec=abs(P1);
lowAmp = find(ampSpec./max(ampSpec) < 0.01);
ampSpec(lowAmp) = 0; 
phaseSpec=atan2(imag(P1),real(P1));

end

