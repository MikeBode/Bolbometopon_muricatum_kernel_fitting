% clear all
% 
% load SDM_1 SDM_1
% load SDM_2 SDM_2
% load SDM_3 SDM_3
% load SDM_4 SDM_4
% load SDM_5 SDM_5
% load SDM_6 SDM_6
% load SDM_7 SDM_7
% load SDM_8 SDM_8
% load SDM_9 SDM_9
% load SDM_10 SDM_10
% load SDM_11 SDM_11
% load SDM_12 SDM_12
% load SDM_13 SDM_13
% load SDM_14 SDM_14
% load SDM_15 SDM_15
% 
% NDistanceMatrix = [SDM_1 ...
%                   SDM_2 ...
%                   SDM_3 ...
%                   SDM_4 ...
%                   SDM_5 ...
%                   SDM_6 ...
%                   SDM_7 ...
%                   SDM_8 ...
%                   SDM_9 ...
%                   SDM_10 ...
%                   SDM_11 ...
%                   SDM_12 ...
%                   SDM_13 ...
%                   SDM_14 ...
%                   SDM_15];
% 

clear all
load StoredDistanceMatrix.mat

C = 1;
STEP = 250;
for c = 1:62
    eval(['SDM_' num2str(c) ' = DistanceMatrix(:,C:C+STEP-1);']);
    eval(['save SDM_' num2str(c) ' SDM_' num2str(c)])
    C = C + STEP;
end
eval(['SDM_' num2str(c) ' = DistanceMatrix(:,(C-STEP):end);']);
eval(['save SDM_' num2str(c) ' SDM_' num2str(c)])


