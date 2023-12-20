% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

%%

function  [resGrdOut,indGrd] = gridTestVars(varargin)
% varargout =
% Outputs N-dimensional grids of test vectors and results of dimension N.
% That last input of varargin should be the vector of dependent results
% (i.e., for power as a function of PID controller gains, varargin =
% [Kp,Ki,Kd,Power].)

% All input vectors must be the same length but do not need to contain the
% same numbers of unique elements

%% check vector lengths
for k = 1:nargin
    L(k)=length(varargin{k});
end
minL=min(L);
maxL=max(L); 
if minL-maxL ~= 0
    error('All vectors must be the same length you scrub')
end

%% find unique elems of input vectors 

for k=1:nargin-1
    unqInd{k}(:,1)= unique(varargin{k});
end

%% create meshing command
outputStr='';
outputStr2='';
inputStr='';
for k=1:nargin-1
    if k == nargin-1
        inputStr = [inputStr  'unqInd{' num2str(k) '}' ];
        outputStr= [outputStr 'X' num2str(k)];
        outputStr2 = [outputStr2 'iSub' num2str(k)];
    else
        inputStr = [inputStr  'unqInd{' num2str(k) '},' ];
        outputStr= [outputStr 'X' num2str(k) ','];
        outputStr2 = [outputStr2 'iSub' num2str(k) ','];
    end
end

%% make meshes of independent variables
eval(['[' outputStr ']= ndgrid(' inputStr ');']);
% this will return gridded independent parameters named X_1 ... X_(nDim)

%% make mesh of dep variable
szMesh=size(X1); % there will always be an X1, and size(X1) == size(XN)
resGrd = zeros(szMesh); % preallocate results grid.
for k=1:length(varargin{1}) % for each point in results vector 
    for kk=1:nargin-1 % find the mapped point grid
        tQ(k,kk)= varargin{kk}(k); 
    end
    for kk=1:nargin-1
        eval(['idx' num2str(kk) '=find(X' num2str(kk) '==tQ(k,kk))'])
    end
    for kk=1:nargin-2
       % if isempty(idxStr)
       %     idxStrN = ['idx' num2str(kk)];
       % elseif k==2
       %     idxStrN = [idxStr, ',idx' num2str(kk)];
       % end
       if kk==1;
            temp = intersect(idx1,idx2)
       elseif kk==2
           eval(['temp2 = intersect(temp,idx' num2str(kk+1),')']);
       else
           eval(['temp3 = intersect(temp2,idx' num2str(kk+1),')']);
           temp2 = temp3
       end
    end
    gridLin=temp2;
    eval(['[' outputStr2 ']= ind2sub(szMesh,gridLin)']);
    eval(['resGrd(' outputStr2 ')=varargin{end}(k)']); 
    clear temp temp2 temp3 idx* iSub*

     % find the corresponding point in results vector and map to results grid 
end

for k=1:nargin-1
    eval(['indGrd{k}=X' num2str(k)]);
end

resGrdOut = resGrd;
indGrdOut = indGrd;

end

