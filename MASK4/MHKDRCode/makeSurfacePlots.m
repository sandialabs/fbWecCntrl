% Copyright 2024 National Technology & Engineering Solutions of Sandia, 
% LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the 
% U.S. Government retains certain rights in this software.
%
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 


function makeSurfacePlots(nDim,pivot,tol,labels,varargin)

% this function creates a contourf plot for N-Dimensional data with each 3D
% plot taken about a pivot point of the remaining dimensions.

%% inputs
% nDim : the number of independent dimensions of the space. For instance,
%   power plotted function of Kp, Ki gains uses nDim = 2. For power P as a
%   function of Kp, Ki, Kd, nDim = 3
% pivot: a point among the independent dimensions about which surfaces will
%   be generated for nDim > 2. For PID example, the pivot will be a Kp,
%   Ki, and Kd point at which the surface will be evaluated (respectively)
%   for P(Ki,Kd), P(Kp,Kd), P(Kp,Ki).
% labels: a string array of length nDim+1 where each entry is an axis label
%   corresponding to a respective varargin (e.g. labels
%   ={'Kp','Ki','Kd','Power'};
% varargin: nDim + 1 separate arguments, where the first nDim arguments are
%   the parameter vectors (e.g., Kp, Ki, Kd) and the additional is the
%   dependent variable vector (e.g., power).

%% outputs
% none explicit, but this will create N figures where N is the number of
% dimesions (except if N = 2, then just 1 figure.

%% check inputs
if nDim ~= nargin-5
    error('nDim does not agree with number of arguments supplied.')
elseif length(pivot) ~=nargin-5
    error('length(pivot) does not agree with number of arguments supplied')
end

%% mesh dependent data
[resGrd,indGrd] = gridTestVars(varargin{:});

%% make surfaces about each pivot
P = nchoosek([1:nDim],2); % permute the number of dimensions. this always takes the plot at the pivot point of indices >2
idxStr = '';
outputStr = '';
szMesh = size(indGrd{1});
for k=1:nDim;
    if k==nDim
        idxStr = [idxStr ':'];
        outputStr = [outputStr 'iSub' num2str(k)];
    else
        idxStr= [idxStr ':,'];
        outputStr = [outputStr 'iSub' num2str(k) ','];
    end
end

idxStrDef = idxStr;
for k=1:length(P(:,1))
    idxStr = idxStrDef;
    exc = 1:nDim;
    exc(P(k,:))=[];
    tLabs = labels(P(k,:)); % labels of plotted values
    tPivot = pivot(exc); % pivot about excluded points
    tVarg{1} = indGrd{P(k,1)};
    tVarg{2} = indGrd{P(k,2)};

    for kk=1:length(exc); % first two dimensions are going to be plotted
        tPivotGrd(kk) = find(abs(indGrd{exc(kk)}(:)-tPivot(kk)) < tol,1,'first'); % this will be index along kk-th dimension containing pivot points
        eval(['[' outputStr ']= ind2sub(szMesh,tPivotGrd(kk))']);
        idxStr = insertAfter(idxStr,2*exc(kk)-1,eval(['num2str(iSub' num2str(exc(kk)) ')']));
        idxStr(2*exc(kk)-1) = '';
        %end
    end

    clear tVarg
    % input string
    tVarg{1} = eval(['squeeze(indGrd{P(k,1)}(' idxStr '))']);
    tVarg{2} = eval(['squeeze(indGrd{P(k,2)}(' idxStr '))']);
    tResGrd = eval(['squeeze(resGrd(' idxStr '))']);

    figure; clf;
    contourf(tVarg{1},tVarg{2},tResGrd)%,11,'LevelList',[30:(max(max(tResGrd))-30)/10:max(max(tResGrd))]);
    colormap(gca,'parula');
    hold on; grid on;
    xlabel(tLabs{1})
    ylabel(tLabs{2})
    c= colorbar;
    c.Label.String = labels{end};
   % clim([30 max(max(tResGrd))]);
   
    clear tVarg idxStr tResGrd exc tLabs pLast pNew tPivotGrd
end
end