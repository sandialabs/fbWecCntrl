% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

%%

function varargout = importNSWCCD(testNum)
% importNSWCCD(testNum)
%
% data = importNSWCCD(testNum)
% Returns structure with fields for each NSWCCD DAQ
%
% [BADAQ,BPDAQ,SAADAQ] = importNSWCCD(testNum)
% Returns seperate table objects.
%
% Call with NSWCCD "Wave Series" number (column A in NSWCCD test log).
%--------------------------------------------------------------------------


daqName = {'BPDAQ','BADAQ','SAADAQ'};

formatSpec{1} = '%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
formatSpec{2} = '%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
formatSpec{3} = '%*s%f%f%f%f%f%f%f%[^\n\r]';

% variable names for each DAQ 
varNames{1} = {'UTCTime','Time','BRIDGEPROBE03','BRIDGEPROBE04','BRIDGEPROBE05','BRIDGEPROBE06','BRIDGEPROBE08','DIRECTIONALARRAY01','DIRECTIONALARRAY02','DIRECTIONALARRAY03','DIRECTIONALARRAY04','DIRECTIONALARRAY05','DIRECTIONALARRAY06','DIRECTIONALARRAY07','DIRECTIONALARRAY08','DIRECTIONALARRAY09','DIRECTIONALARRAY10','DIRECTIONALARRAY11','DIRECTIONALARRAY12','WAVEMAKERPULSE','SINEWAVE','STARTCOLLECTTTL','BUOY01','BUOY02','BUOY03','BUOY04','BUOY05'};
varNames{2} = {'UTCTime','Time','OSSI01','OSSI02','OSSI03','OSSI04','OSSI05','OSSI06','OSSI07','OSSI08','OSSI09','OSSI10','OSSI11','OSSI12','StringPot'};
varNames{3} = {'UTCTime','Time','SAA01','SAA02','SAA03','SAA04','SAA05'};

% mark the sensors that are wave probes and should be scaled (inches to m)
% and run through spectral analysis. Note that these reference the
% timetable inds (1 off from varNames)
isWaveProbe{1} = [3:19,23:27];
isWaveProbe{2} = 3:14;
isWaveProbe{3} = 3:7;

% loop through the three DAQs (BPDAQ, BADAQ, SAADAQ)
for ii = 1:length(daqName)
    try
        [mat{ii},filename] = readTwe(daqName{ii},testNum,formatSpec{ii},varNames{ii},isWaveProbe{ii});
        data.(daqName{ii}) = mat{ii};
        data.src{ii} = filename;
        
        % find wave maker start time
        if ii == 1 && strcmp(daqName{ii},'BPDAQ')
            data.BPDAQ = findWaveMakerPulse(data.BPDAQ);
            data.waveMakerStartTime = data.BPDAQ.Time(find(data.BPDAQ.waveMakerOn,1));
        end
        
%         % calculate spectral properties
%         [spectra{ii},chars{ii}] = getSpecChars(data.(daqName{ii}),data.BPDAQ,isWaveProbe{ii});
        
    catch ME
        warning(getReport(ME))
        warning('%s not found',daqName{ii});
        data.(daqName{ii}) = [];
    end
end

% data.chars = vertcat(chars{:});
% data.spectra = mergeStructs(spectra);

% Formulate output
if nargout == 1 || nargout == 0
    
    % structure with Timetables
    varargout = {data};
elseif nargout == 3
    
    % matrices
    varargout = mat;
end

end


function BPDAQ = findWaveMakerPulse(BPDAQ)
% BPDAQ = findWaveMakerPulse(BPDAQ)

waveMakerOn = 1*(BPDAQ.WAVEMAKERPULSE>1);
BPDAQ = addvars(BPDAQ,waveMakerOn);
end


function [spectra, chars] = getSpecChars(ttab,BPDAQ,isWaveProbe)
% Use WAFO to compute spectra and spectral characteristics 

% find when wave maker is on
tmp = interp1(BPDAQ.Time,BPDAQ.waveMakerOn,ttab.Time,'nearest',0);
inds = find(tmp);

jj = 0;
for ii = isWaveProbe-1
    jj = jj + 1;
    S = dat2spec([ttab.Time(inds),ttab{inds,ii}],400);
    [chartmp(:,jj),~,chtext] = spec2char(S,1:15);
    spectra.(ttab.Properties.VariableNames{ii}) = S;
end

chars = array2table(chartmp');
chars.Properties.RowNames = ttab.Properties.VariableNames(isWaveProbe-1);
chars.Properties.VariableNames = chtext;
chars.Properties.VariableUnits = {'m','s','s','s','s','s','','','','','s','','','',''};

end


function fname = getfname(testNum,daqName)
% fname = getfname(testNum,daqName)
%
% Construct filenames based on test number and DAQ name

fname = sprintf('AWDC-%s-s%05d.twe',daqName,testNum);
if ~exist(fname,'file')
    error('No such files in path: %s',fname)
end
end


function [ttab,filename] = readTwe(daqName,testNum,formatSpec,varNames,isWaveProbe)
% [ttab,filename] = readTwe(daqName,testNum,formatSpec,varNames,isWaveProbe)
%
% Parse TWE files and create a Timetable structure

delimiter = '\t';
startRow = 20;
endRow = inf;
filename = getfname(testNum,daqName);
fileID = fopen(filename,'r');
fprintf('Reading from %s\n',filename)
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

tab = table(dataArray{1:end-1}, 'VariableNames',varNames);

frewind(fileID)
D = textscan(fileID,'%s\t%f',1,'headerlines',1);
dd = D{2};
T = textscan(fileID,'%s\t%s',1,'headerlines',1);
tt = T{2}{1};
stime = datetime([num2str(dd),'-',tt],'InputFormat','yyyyMMdd-HH:mm:ss');

if tab.UTCTime(1) == tab.UTCTime(2)
    tab.UTCTime = stime + seconds(tab.Time);
else
    tab.UTCTime = datetime(num2str(tab.UTCTime),'InputFormat','HHmmss.SSSSSS');
end

ttab = table2timetable(tab);


% scale from inches to meters
for ii = isWaveProbe - 1
    ttab{:,ii} = ttab{:,ii}*0.0254;
end

fclose(fileID);
end

function [s] = mergeStructs(sca)

for ii = 1:length(sca)
    fn = fields(sca{ii});
    for jj = 1:length(fn)
        s.(fn{jj}) = sca{ii}.(fn{jj});
    end
end

end
