function OutputName = saveImage(Mat, V, OutputName);
%
% Syntax :
%     OutputName = saveImage(Mat,V,OutputName);
%
% This function saves a 4D matrix in as a nifti file.
%
% Input Parameters:
%        Mat                    : 4D Matrix
%        V                      : Volume estructure
%        OutputName             : Output Filename
%
% Output Parameters:
%        OutputName             : Output Filename
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez and Javier Santoja
% LIM, HUGGM
% November 13th 2014
% Version $1.0

for j = 1:size(Mat,4)
    Vol(j) = V;
    Vol(j).fname = OutputName;
    Vol(j).dt = [spm_type('float32') 0];
    Vol(j).descrip = sprintf('%s','');
    Vol(j).n = j;
end
Vol =spm_create_vol(Vol);
for j =1:size(Mat,4);
    spm_write_vol(Vol(j),squeeze(Mat(:,:,:,j)));
end
return