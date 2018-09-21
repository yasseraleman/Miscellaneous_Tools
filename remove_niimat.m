function remove_niimat(varargin);
%
% Syntax :
%      remove_niimat(niiFname);
%
% This script removes the mat file created  during the nifti saving
% process.
%
% Input Parameters:
%        niiFname                  : Image filename or Image volume (see spm_vol);
%
% Output Parameters:
%
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% Cuban Neuroscience Center
% February 21th 2006
% Version $1.0

niiFname = varargin{1};
[pth, nm, ext] = fileparts(niiFname);

if exist([pth filesep nm filesep '.mat']);
    delete([pth filesep nm filesep '.mat']);
end
return;