function outPaths = Path_Extraction(InputFolder);
%
% Syntax :
% outPaths = Path_Extraction(InputFolder);
%
% This function extracts all the "-study-" paths from a specified input folder.
%
% Input Parameters:
%      InputFolder   : Input Folder
%
%
% Output Parameter:
%     outPaths      : Extracted Paths
%
% See also: Surf_Color
%__________________________________________________
% Authors: Yasser Aleman Gomez
%LIM
% February 1th 2013
% Version $1.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TempDir = genpath(InputFolder);
if isunix
    Pthst = strread(TempDir,'%s','delimiter',':'); Pthst = char(Pthst);
else
    Pthst = strread(TempDir,'%s','delimiter',';'); Pthst = char(Pthst);
end
outPaths = '';
cont = 0;
indp = '';
for i = 1:size(Pthst,1)
    ind=strfind(lower(Pthst(i,:)),'-study-');
    if ~isempty(ind)
        cont = cont+1;
        indp(cont) = i;

    end
end
if ~isempty(indp)
    outPaths = Pthst(indp,:);
end
for i =1:size(outPaths,1);ind = strfind(outPaths(i,:),filesep);a(i) = length(ind);end
ind = find(a==max(a));outPaths = outPaths(ind,:);
return