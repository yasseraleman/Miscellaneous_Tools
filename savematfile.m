function [Outfile] = savematfile(Outfile, matfile);
% Syntax :
% [Outfile] = savematfile(Outfile, matfile);
%
% This function saves matlab variables with a specified filename.
%
% Input Parameters:
%   Outfile     :  Name for the output matlab file.
%   matfile     :  Matlab variable name.
%
% Output Parameters:
%   Outfile     :  Name for the output matlab file.
%
% Related references:
% 
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

%=========================Main program====================================%
[npth,nm,ext] = fileparts(Outfile);

if isstruct(matfile)&isfield(matfile,'SurfData')
    Surf = matfile;
    VarName = ' Surf';
elseif isstruct(matfile)&isfield(matfile,'Vnorm')
    ResultsVN = matfile;
    VarName = ' ResultsVN';
elseif isstruct(matfile)&~isfield(matfile,'Area')
    ResultsTV = matfile;
    VarName = ' ResultsTV';
elseif isstruct(matfile)&isfield(matfile,'Area')
    ResultsVA = matfile;
    VarName = ' ResultsVA';
else
    data =whos('matfile');
    VarName =[' ' data.name];
end
inds = isspace(nm);
if ~isempty(inds)
    nm(inds)= '_';
end
if sum(isspace(Outfile))~=0
    if ~isunix
        eval(['save ' [npth(1:3) nm '.mat'] VarName]);
        movefile([npth(1:3) nm '.mat'],npth);
        Outfile = [npth filesep nm '.mat'];
    elseif isunix
        [stat,cdir] = system('echo $HOME');
        eval(['save ' [cdir nm '.mat'] VarName]);
        movefile([cdir nm '.mat'],npth);
        Outfile = [npth filesep nm '.mat'];
    end
else
    eval(['save ' deblank(Outfile) VarName]);
end
%========================End of main program==============================%
return;