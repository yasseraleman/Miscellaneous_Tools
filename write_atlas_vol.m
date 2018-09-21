function V = write_atlas_vol(V,Y);
% Syntax :
% V = write_atlas_vol(V,Y);
%
% This function is based on spm_write_vol function developed by PhD. John Ashburner(FIL,UCL). 
% It has been adapted to keep same scales factor and offsets of original
% image V.
%
% Input Parameters:
%   V     :  A structure containing image volume information (see spm_vol).
%   Y     :  A one, two or three dimensional matrix containing the image
% voxels
%
%
% Output Parameters:
% V (output) - data structure after modification for writing.
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

if ndims(Y)>3, error('Can only handle a maximum of 3 dimensions.'), end

if ~isfield(V,'pinfo'), V.pinfo = [1,0,0]'; end

dim = [size(Y) 1 1 1];
if ~all(dim(1:3) == V.dim(1:3)) | (size(V.pinfo,2)~=1 & size(V.pinfo,2)~=dim(3)),
	error('Incompatible dimensions.');
end

% Set scalefactors and offsets
%-----------------------------------------------------------------------
V.pinfo(1,:) = 1;
V.pinfo(2,:) = 0;

%-Create and write image
%-----------------------------------------------------------------------
V = spm_create_vol(V);
for p=1:V.dim(3),
	V = spm_write_plane(V,Y(:,:,p),p);
end;
fclose all;