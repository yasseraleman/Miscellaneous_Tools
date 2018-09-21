function varargout = crop_Images(imageFilename, varargin);
%
% Syntax :
%      croppedImage = crop_Images(imageFilename, borders, boolsave);
%
% This script crops an image to remove empty voxels in its borders. Cropping image
% can reduce make computation time, and memory consumption.
%
% Input Parameters:
%        imageFilename             : Image filename or Image volume (see spm_vol);
%        borders                   : Crop Limits
%        boolsave                  : Boolean variable to save or not the
%                                    image.
%
% Output Parameters:
%        croppedImage             : Cropped Image Filename.
%
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% Cuban Neuroscience Center
% February 21th 2006
% Version $1.0

%% ====================== Checking input parameters ===================== %
if nargin<1 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    % Reading ODF File or Variable
    temp = whos('imageFilename');
    switch temp.class
        case 'char'
            if ~exist(imageFilename,'file')
                error('The image file does not exist');
                return
            else
                [pth,nm,ext] = fileparts(imageFilename);
                switch deblank(ext)
                    case '.gz'
                        boolzipImage = 1;
                        tempName = unzip_nifti(imageFilename);
                    case '.nii'
                        boolzipImage = 0;
                        tempName = imageFilename;
                    otherwise
                        error('Unrecognized image format.');
                        return;
                end
                V = spm_vol(tempName);
                
                param.Xlim = [1 V(1).dim(1)];
                param.Ylim = [1 V(1).dim(2)];
                param.Zlim = [1 V(1).dim(3)];
                param.Tlim = [1 length(V)];

                voxsize = sqrt(sum(V(1).mat(1:3,1:3).^2));
            end
        case 'double'
            Ptemp = imageFilename;
            param.Xlim = [1 size(I,1)];
            param.Ylim = [1 size(I,2)];
            param.Zlim = [1 size(I,3)];
            param.Tlim = [1 size(I,4)];
        case 'struct'
            V = imageFilename;
            tempName = V(1).fname;
            param.Xlim = [1 V(1).dim(1)];
            param.Ylim = [1 V(1).dim(2)];
            param.Zlim = [1 V(1).dim(3)];
            param.Tlim = [1 length(V)];
    end

    % Parameters
    boolSave = 1; % Boolean Variable to save or not the cropped image.
end


% deal with the input arguments
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'xLims'
                    % Limits in X axis (Left-Right axis)
                    param.Xlim=varargin{2};
                case 'yLims'
                    % Limits in X axis (Anterior-Posterior axis)
                    param.Ylim=varargin{2};
                case 'zLims'
                    % Limits in Z axis (Inferior-Superior axis)
                    param.Zlim=varargin{2};
                case 'tLims'
                    param.Tlim=varargin{2};
                case 'boolSave' % Boolean Variable to save or not the cropped image.
                    boolSave=varargin{2};
                    if isnumeric(boolSave)
                        boolSave = logical(boolSave);
                    else
                        error('Wrong boolean variable');
                        return
                    end
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end

if nargout > 1
    errordlg('To Many Output Parameters');
    return;
end

borders = [param.Xlim;param.Ylim;param.Zlim;param.Tlim];

%% ========================= End of Checking Inputs ==================== %%

try
    
    %% ======================== Main Program ================================ %
    cont = 0;
    cropDim = borders(:,2) - borders(:,1) +1;
    for j = borders(4,1):borders(4,2)
        cont = cont  + 1;
        if exist('V','var')
            Ptemp = spm_read_vols(V(j));
        end
        P(:,:,:,cont) = Ptemp(borders(1,1):borders(1,2),borders(2,1):borders(2,2),borders(3,1):borders(3,2));
        cropDim = size(P);
        if length(cropDim) == 2
            cropDim(3) = 1;
        end
        %% ======================== Saving Cropped Image ==================== %
        
        if boolSave&exist('V','var')
            Volt(cont) = V(1);
            [pth,nm,ext] = fileparts(V(1).fname);
            Volt(cont).fname = fullfile(pth,[nm '_crop.nii']);
            Volt(cont).dt = [spm_type('float32') 0];
            Volt(cont).descrip = sprintf('%s','Cropped Image');
            Volt(cont).n = cont;
            Volt(cont).dim(1:3) = [cropDim(1) cropDim(2) cropDim(3)];
            vec = [borders(1,1)-1 borders(2,1)-1 borders(3,1)-1];
            vec = Volt(cont).mat(1:3,1:3)*vec';
            Volt(cont).mat(1:3,4) = Volt(cont).mat(1:3,4) + vec;
        end
        
    end
    if boolSave&exist('Volt','var')
        Volt =spm_create_vol(Volt);
        for j = 1:size(P,4)
            spm_write_vol(Volt(j),P(:,:,:,j));
        end
        remove_niimat(Volt(1).fname);
        if boolzipImage
            outFilenamezip = gzip(Volt(1).fname);
            delete(Volt(1).fname);
            varargout{1} = outFilenamezip;
            gzip(V(1).fname);
            delete(V(1).fname);
        else
            varargout{1} = Volt(1).fname;
        end
    else
        varargout{1} = P;
    end
catch
    switch ext
        case '.gz'
            gzip(tempName);
    end
    varargout{1} = imageFilename;
end

%% ======================= End of Main Program ========================== %
return;