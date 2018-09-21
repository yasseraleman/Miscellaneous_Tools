function varargout = Create_4D_SubVolumes(varargin);
%
% Syntax :
% OutFilename = Create_4D_SubVolumes(DWIImage, Rindex, OutFilename);
%
% This script creates a new 4D volume just using the specified volumes in
% Rindex
%
% Input Parameters:
%   DWIImage          : 4D Volume
%   Rindex            : Volumes indexes to create the new 4D Volume
%
% Output Parameters:
%   OutFilename       : New 4D image
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 22th 2017
% Version $1.0

%% ============================ Checking input parameters ============================== %
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
else
    DWIImage = varargin{1};
    Rindex = varargin{2};
    [outTempDir,outTempFile] =  fileparts(DWIImage);
    
    OutFilename = [outTempDir filesep outTempFile '_subVol.nii'];
end

if nargin == 3
    OutFilename = varargin{3};
end

if nargin > 3
    error('To many input arguments');
    return;
end
%% ====================== End of Checking input parameters ============================= %

%% ========================== Creating New 4D volumes ================================== %
if ischar(Rindex)
    eval(['ind2load = ' Rindex(:)']);
else
    ind2load = Rindex(:)';
end

Vdiff = spm_vol(DWIImage);
V = Vdiff(1);
for j = 1:length(Rindex)

    Vol(j) = V;
    Vol(j).fname = OutFilename;
    Vol(j).descrip = sprintf('%s','4D Subvolume');
    Vol(j).n = [j 1];
end
Vol = spm_create_vol(Vol);


Ns = length(ind2load);
for i = 1:Ns
%     disp(['Writing new volume: Number ' num2str(i) ' of ' num2str(Ns)]);
    I = spm_read_vols(Vdiff(ind2load(i)));
    spm_write_vol(Vol(i),I);
end

%% ========================== End of Creating New 4D volumes =========================== %
% Outputs
varargout{1} = OutFilename;
return;