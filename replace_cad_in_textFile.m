function varargout = replace_cad_in_textFile(varargin);
%
% Syntax :
%  [groupsIds, groupsNames, allsulciIds, labelsId] = Detecting_BrainVisa_Groups(ChangeTxtFile, hemi);
%
% This script creates a stats table using the BrainVisa outputs.
%
% Input Parameters:
%         inFile                    : Input text filename.
%         outFile                   : Output text filename.
%         inCad                     : String that will be replaced.
%         outCad                    : New string.
%
%
% Output Parameters:
%       outFile                     : Output text filename.
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin <4
    error('Four inputs are mandatory');
    return;
end

inFile = varargin{1};
outFile = varargin{2};
inCad = varargin{3};
outCad = varargin{4};
%% ======================= End of Input parameters  ======================%



fin = fopen(inFile);
fout = fopen(outFile,'wt');

while ~feof(fin)
   s = fgetl(fin);
   s = strrep(s, inCad, outCad);
   fprintf(fout,'%s\n',s);
end

fclose(fin);
fclose(fout);
varargout{1} = outCad;
return;