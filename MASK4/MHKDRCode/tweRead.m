% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 


function tweMat = tweRead(tweFileName);

% reads TWE file and builds *mat structure of outputs
% Input: a string containing the full file name (including path) to the 
%   *.twe file.'
% Output: a data structure containing the basin wave sensors

rawT = readtable(tweFileName,'NumHeaderLines',5,'FileType','text')
dataStartRow = find(strcmp(rawT.CHANNELNAME_,'DATA:'));
rawT = rawT(dataStartRow:end,2:end); % first column of data will be NaNs
rawT{:,[4:22]} = detrend(rawT{:,[4:22]},'constant') * 0.0254;

%% convert time stamps
fileID = fopen(tweFileName);
D = textscan(fileID,'%s\t%f',1,'headerlines',1);
dd = D{2};
T = textscan(fileID,'%s\t%s',1,'headerlines',1);
tt = T{2}{1};
stime = datetime([num2str(dd),'-',tt],'InputFormat','yyyyMMdd-HH:mm:ss');
rawT.UTCTime = datetime(string(stime,'yyyy-MM-dd-') ...
    + num2str(rawT.UTCTime),'InputFormat','yyyy-MM-dd-HHmmss.SSSSSS');
%rawT.TAITime = datetime(rawT.TAITime/1e9,'convertfrom','posix');
tweMat = table2struct(rawT,'ToScalar',true);
fclose(fileID);

end