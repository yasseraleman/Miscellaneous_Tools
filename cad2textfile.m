function varargout = cad2textfile(varargin);
%
% Syntax :
%   textFile = cad2textfile(cad2print, textFile);
%
% This script writes lines included from the variable cad2print into
% a text file.
%
% Input Parameters:
%         cad2print        : Lines to be printed.
%         textFile         : Output text file.
%
% Output Parameters:
%         functNames       :  Functions filenames
%
%
% Related references:
%
%
% See also: 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 14th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%


%% ========================== Checking Inputs ========================== %%
if nargin <2
    error('Two inputs are mandatory');
    return;
end
cad2print = varargin{1};
textFile = varargin{2};
%% ====================== End of Checking Inputs ======================= %%

%% ========================== Main Program ============================= %%
fid = fopen(textFile,'wt');
for i = 1:size(cad2print,1)
    fprintf(fid,'%s\n',cad2print(i,:));
    
end
fclose(fid);
%% ======================= End of Main Program ========================= %%

% ---- Outputs
varargout{1} = textFile;
return