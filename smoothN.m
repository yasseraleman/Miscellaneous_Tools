function g = smoothN(f,scales,ox,oy,oz,equal_cont)

%SMOOTHN smoothes a multivalued (scalar, vector, tensor or any higher order)
% image, or its derivative, with a gaussian kernel.
%
% USAGE: I = SMOOTHN(IMAGE,SCALES,DX,DY,DZ,EQ_CONT)
%
% IMAGE: image to be smooth. IMAGE can be an array from the workspace, an
% SPM volume or a string with the filename.
%
% SCALES: widths of the gaussian kernel. SCALES is a vector of length 3 to
% specify width of the gaussian kernel in space.
%
% DX, DY, DZ: options to specify derivative cases. they can be 0, 1 or 2
% depending on the order of the derivative.
%
% EQ_CONT: special option to set the gaussian kernel to the mode of equal
% contribution. SCALES, DX, DY and DZ options will be disabled
% and the shape of the kernel along each dimension takes the form:
% [.05 .25 .4 .25 .05].
%
% I: array of the smoothed image.
%
% See also: GETWARPS, DOWARP
%
% Author: Pedro Antonio Valdes-Hernandez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 18 2005
% Version $1.0

%=========================Main program====================================%
if isstruct(f)
    f = spm_read_vols(f);
elseif ischar(f)
    f = spm_read_vols(spm_vol(f));
end

S = size(f);
g = zeros(S);

if scales == [0 0 0];
    g = f;
else
    if equal_cont
        Gsx = [.05 .25 .4 .25 .05];
        Gsy = Gsx; Gsz = Gsx; Kz = 2;
    else
        % Calculate sample points
        Kx = ceil(scales(1));
        x = -Kx:Kx;
        Ky = ceil(scales(2));
        y = -Ky:Ky;
        Kz = ceil(scales(3));
        z = -Kz:Kz;

        % Sample Gaussian function and normalize
        Gsx = exp(-x.^2/(2*scales(1)^2));
        Gsx = Gsx / sum(Gsx);
        Gsy = exp(-y.^2/(2*scales(2)^2));
        Gsy = Gsy / sum(Gsy);
        Gsz = exp(-z.^2/(2*scales(3)^2));
        Gsz = Gsz / sum(Gsz);

        % Calculate the derivates in x and y directions
        Gsx = derivativeN(ox,x,Gsx,scales(1));
        Gsy = derivativeN(oy,y,Gsy,scales(2));
        Gsz = derivativeN(oz,z,Gsz,scales(3));
    end

    % Do the convolutions
    for ind = 1:prod(S(3:end))
        g(:,:,ind) = conv2(Gsx,Gsy,f(:,:,ind),'same');
    end
    g(:,:,end + (2 * Kz),:) = 0;
    g = filter(Gsz,1,g,[],3);
    g(:,:,[1:Kz end-Kz + 1:end],:) = [];
end
%========================End of main program==============================%
return

%========================Internal Functions===============================%
function r = derivativeN(order,x,Gs,scale)

switch order
    case 0
        r = Gs;
    case 1
        r = -(x/(scale^2)).*Gs;
    case 2
        r = (x.^2 - scale^2)/(scale^4).*Gs;
    otherwise
        error('only derivatives up to second order are supported');
end
return;