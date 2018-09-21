function  varargout = Remove_Inserted_Spaces(varargin);
%
% Syntax :
%    outNames = Remove_Inserted_Spaces(inNames);
%
% This function removes white space inserted in different input names
%
% Input Parameters:
%     inNames                 : Input Names.
%
% Output Parameters:
%     outNames                : Output Names.
%
%
% Related references:
%
% See also: 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% August 8th 2012
% Version $1.0

%% ================== Checking Input parameters ========================= %
if nargin <1
    error('Please enter a correct number of inputs');
    return;
end
inNames = varargin{1};
if ischar(inNames)
    tempNames = inNames;
elseif iscell(inNames)
    tempNames = char(inNames);
elseif iscell(inNames)|ischar(inNames)
    error('The input Names are not cell or char arrays');
    return;
end
outNames = '';
for i = 1:size(tempNames,1);
    temp = deblank(tempNames(i,:));
    ind = isspace(temp);
    temp(ind) = [];
    outNames = strvcat(outNames, deblank(temp));
end
if iscell(inNames)
    outNames = cellstr(outNames);
end
varargout{1} = outNames;
return;