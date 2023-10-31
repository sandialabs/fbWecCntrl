function [t,f,OUT] = genMultiSine_NInput(fmin,fmax,T,options)
% genMultiSine   generate multisine signal
%
% Generates a set of white or pink multi-sine signals for use in system
% identification (SID) experiments.
%
% A large set (specified by numPhases) of phase realizations is generated
% and the signals are sorted in ascending order of peak-to-peak amplitude;
% signals with peak-to-peak (p2p) larger than 10% of min p2p are discarded.
%
% From the remaining signals, a combination of (NumExp x NumInput) signals
% is sampled. The condition number for inversion of the spectrum at
% each frequency k for each (NumExp x NumInput) group are calculated.
% Lastly, the group with the smallest max conditional number is taken as
% the set of output signals.
%
% Arguments:
%   fmin      minimum frequency [Hz]
%   fmax      maximum frequency [Hz]
%   T         period of signal [s]
%
% optional name-value pairs
%   numPhases   number of phasings to consider (default: 1e2)
%   plotFlag    set to 1 for plots (default: 0)
%   color       'pink' or 'white' (default: 'white')
%   dt          time-step for time series [s] (default: 2e-2)
%   RMS         root-mean-square amplitude (default: 1)
%   k           decay factor for pink noise (default: 0.5)
%   NumExp      number of experiment (default: 1)
%   NumInput    number of input (default: 1)
%   NumRepeat   number of cycle (default: 3)
%
%
% Returns:
%   t         time series vector [s]; [row,col] = [Nt, 1]
%   f         frequency vector [Hz]; [row,col] = [Nf, 1]
%   OUT       NumExp structure:
%                xt        signal time vector [*]; [row,col] = [Nt,NumInput]
%                Xf        signal complex frequency vector [*]; [row,col] = [Nf,NumInput]
%
% Example:
% [t,f,OUT] = genMultiSine_NInput(0.05,2.5,40,...
%    'NumExp',4,'NumRepeat',5,'k',0.75,'color','pink','dt',0.001,...
%    'numPhases',10000,'plotFlag',1);
%
% Authors: Giorgio Bacelli (original), Ryan Coe (minor tweaks),
%          Alicia (edited for N-input system)
%
% -------------------------------------------------------------------------

arguments
    fmin (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
    fmax (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
    T (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
    options.numPhases (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 1e2
    options.plotFlag (1,1) {mustBeNumericOrLogical} = 0
    options.color (1,:) string = 'white'
    options.dt (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 1e-2
    options.RMS (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 1
    options.k (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 0.5
    options.NumExp (1,1) {mustBeFinite,mustBeReal,mustBePositive} = 1
    options.NumInput (1,1) {mustBeFinite,mustBeReal,mustBePositive} = 1
    options.NumRepeat (1,1) {mustBeFinite,mustBeReal,mustBePositive} = 3
end

% create time and frequency vector
df = 1/T;
dt = options.dt;
idx_fmin = round(fmin/df);
idx_fmax = round(fmax/df);
f_vec = (idx_fmin:idx_fmax)'*df;
w_vec = 2*pi*f_vec;
Nf = length(f_vec);

% phase randomized between -pi to pi
rng('shuffle')
ph_mat = -pi + (pi--pi).*rand(Nf,options.numPhases);

% create amplitude vector
switch lower(options.color)
    case 'pink'
        Amp = ones(size(f_vec)) ./ f_vec.^options.k;
    case 'white'
        Amp = ones(size(f_vec));
    otherwise
        error('Unkown noise color (choose ''pink'' or ''white'')')
end

% Create scale amplitude vector for RMS = options.RMS
RMS = sqrt(sum(Amp.^2))/sqrt(2);
GAIN = options.RMS/RMS;
Amp = Amp * GAIN;

% Create exp(iwt) and complex amplitude A*exp(i*theta)
t_vec = 0:dt:T-dt;
Nt = length(t_vec);
exp_mat = exp(1i * (w_vec*t_vec));
Amp_f = Amp.*exp(1i*(ph_mat));

signal = real(Amp_f'*exp_mat)';

% Sort signals in ascending p2p
[p2pval,idx] = sort(max(signal)-min(signal));
signal = signal(:,idx);
Amp_f = Amp_f(:,idx);

% Discard signals with p2p 10% greater than min(p2p)
idx = find(p2pval <= p2pval(1)*1.1);
signal = signal(:,idx);
Amp_f = Amp_f(:,idx);

% if size(Amp_f,2) >= 2*(options.NumInput*options.NumExp)
%     signal = signal(:,1:2*(options.NumInput*options.NumExp));
%     Amp_f = Amp_f(:,1:2*(options.NumInput*options.NumExp));
%     disp('*********** Warning: samples further reduced. ***********')
% end

Nsample = size(Amp_f,2)

% Conditional number
% If a matrix is singular, then its condition number is infinite.

% first check if we have enough samples
if size(Amp_f,2) <= (options.NumInput*options.NumExp)
    disp('*********************************************************')
    disp('Not enough samples,')
    disp('Please re-run for different set of randomly generated phases,')
    disp('OR ')
    disp('increase option.numPhases.')
    disp('*********************************************************')
    error('Not enough samples, please increase option.numPhases!')
end

% Generate list of possible combination of NumInput*NumExp signals
SignalSets = nchoosek((1:size(Amp_f,2)),options.NumInput*options.NumExp);

% Calculate the condition number at each frequency and choose the set with
% smallest max condition number (best U^-1)
currMin = 0;
condIdx = 0;
condVec = zeros(1,Nf); % store all condition number
for ii = 1:size(SignalSets,1)
    for jj = 1:Nf
        testSig = reshape(Amp_f(jj,SignalSets(ii,:)),options.NumExp,options.NumInput);
        condVec(jj) = cond(testSig*ctranspose(testSig));
    end
    if max(condVec) <= currMin
        currMin = max(condVec);
        condIdx = ii;
    end
    if condIdx == 0
        currMin = max(condVec);
        condIdx = ii;
    end
end
SignalSetsUsed = reshape(SignalSets(condIdx,:),options.NumExp,options.NumInput);

%%
% New time vector
t = (0:dt:options.NumRepeat*Nt*dt-dt)';
f = f_vec;

OUT =  struct('xt',{},'Xf',{});

for ii = 1:options.NumExp
    OUT(ii).xt = repmat(signal(:,SignalSetsUsed(ii,:)),options.NumRepeat,1);
    OUT(ii).Xf = Amp_f(:,SignalSetsUsed(ii,:));
end

if options.plotFlag

    figure
    stem(f_vec, Amp, '.')
    xlabel('Frequency [Hz]','interpreter','latex');
    ylabel('Amplitude','interpreter','latex');
    title(sprintf('%s noise, $f \\in [%.2g,%.2g]$, RMS = %g',...
        options.color,fmin,fmax, options.RMS),'interpreter','latex');
    grid on


    for ii = 1:options.NumExp
        figure
        for jj = 1:options.NumInput
            subplot(options.NumInput,1,jj)
            plot(t, OUT(ii).xt(:,jj))
            hold on
            plot(t_vec,real(OUT(ii).Xf(:,jj)'*exp_mat)','--')
            legend('OUT','calculation')
            xlabel('Time [s]','interpreter','latex');
            ylabel('Ampiltude','interpreter','latex');
            title(['Signal: ' num2str(jj)])
            grid on
        end
        sgtitle(['Experiment: ', num2str(ii)],'interpreter','latex');
    end
end

end


