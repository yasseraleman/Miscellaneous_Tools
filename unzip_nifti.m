function varargout = unzip_nifti(varargin);
%
% Syntax :
%      zipFilename = unzip_nifti(zipFilename, outDir);
%
% This script decompress nifty files.
%
% Input Parameters:
%        zipFilename               : Image filename.
%        outDir                    : Output directory.
%
%
% Output Parameters:
%        uzipFilename               : Unzipped Image.
%
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% Cuban Neuroscience Center
% March 22th 2007
% Version $1.0

%% ============================= Checking Inputs ======================= %%

if nargin == 1
    zipFilename = varargin{1};
    [outDir, nm, ext] = fileparts(zipFilename);
    boolsave = 0;
elseif nargin == 2
    zipFilename = varargin{1};
    outDir      = varargin{2};
    boolsave = 0;
elseif nargin == 3
    zipFilename = varargin{1};
    outDir       = varargin{2};
    boolsave     = varargin{3};
else
    error('Please check input parameters');
    return;
end
if ~exist(zipFilename,'file')
    error('The specified zip file does not exist');
    return;
else
   [pth, nm, ext] = fileparts(zipFilename);
end
if isnumeric(boolsave)
    boolsave = logical(boolsave);
else
    warning('The boolsave variable must be numeric. The original file will be deleted: boolsave = 0');
    boolsave = 0;
end
%% ======================= End of Checking Inputs ====================== %%

%% ============================ Main Program =========================== %%
switch ext
    case '.gz'
        tempName = [outDir filesep nm];
        [pth, nm, ext] = fileparts(tempName);
        indnii = strfind(nm,'.nii');
        if ~isempty(indnii)
            if indnii == length(nm)-3;
                nm = nm(1:indnii-1);
            end
        end
        gunzip(zipFilename, outDir)
        uzipFilename =    [outDir filesep nm '.nii'];
        
        if ~boolsave
            delete(zipFilename);
        end
        
        %  Output
        varargout{1} = uzipFilename;
    case '.nii'
        disp('This file is not a zipped file.')
        uzipFilename = zipFilename;
        
        %  Output
        varargout{1} = uzipFilename;
end
%% ====================== End of Main Program ========================== %%
return