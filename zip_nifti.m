function varargout = zip_nifti(varargin);
%
% Syntax :
%      zipFilename = zip_nifti(uzipFilename, outDir, boolsave);
%
% This script decompress nifty files.
%
% Input Parameters:
%        uzipFilename              : Image filename.
%        outDir                    : Output directory.
%        boolsave                  : Boolean variable to remove unzipped
%                                    image (1: Keep unzipped file,
%                                    0(default) otherwise).
%
% Output Parameters:
%        zipFilename               : Zipped Image.
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
    uzipFilename = varargin{1};
    [outDir, nm, ext] = fileparts(uzipFilename);
    boolsave = 0;
elseif nargin == 2
    uzipFilename = varargin{1};
    outDir       = varargin{2};
    boolsave = 0;
elseif nargin == 3
    uzipFilename = varargin{1};
    outDir       = varargin{2};
    boolsave     = varargin{3};
else
    error('Please check input parameters');
    return;
end
if ~exist(uzipFilename,'file')
    error('The specified zip file does not exist');
    return;
else
    [pth, nm, ext] = fileparts(uzipFilename);
end
if isnumeric(boolsave)
    boolsave = logical(boolsave);
else
    warning('The boolsave variable must be numeric. The original file will be deleted: boolsave = 0');
    boolsave = 0;
end
%% ======================= End of Checking Inputs ====================== %%

%% ============================ Main Program =========================== %%
switch deblank(ext)
    case '.nii'
        gzip(uzipFilename, outDir);
        uzipFilename =    [outDir filesep nm '.nii'];
        
        %  Output
        varargout{1} = uzipFilename;
        if ~boolsave
            delete(uzipFilename);
        end
    case '.gz'
        disp('This file was previously zipped.')
        zipFilename = uzipFilename;
        
        %  Output
        varargout{1} = zipFilename;
end

%% ====================== End of Main Program ========================== %%
return