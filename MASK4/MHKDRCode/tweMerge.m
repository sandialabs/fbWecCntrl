
% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 


function waveBotFull = tweMerge(tweSync,waveBotMat);

 % inputs
%   tweSync: the data structure of syncd twe file, output of tweSync.
%   waveBotMat: the data structure of the waveBot, output of loadMask4
% outputs: 
%   waveBotFull: waveBotMat with an additional field with wave data

waveBotMat.waveData= tweSync;
waveBotFull = waveBotMat; 

end