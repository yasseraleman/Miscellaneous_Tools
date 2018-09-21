function Build_4DVolumes(InputFolders);
%
% Syntax :
%  Build_4DVolumes(InputFolders);
%
% This function Converts Niftis inside the specified folders list to 4D
% volumes.
%
% Input Parameters:
%   InputFolders  : Folders to construct 4D Nifti files
%
% Output Parameters:
%
% See also: Separate_Dicoms NiftiConv_Dicoms
%__________________________________________________
% Authors:  Yasser Aleman Gomez 
% Laboratorio de Imagenes Medicas.HGGM
% March 3rd 2012
% Version $1.0

if size(InputFolders,1)==1
    TempDir = genpath(InputFolders);
    if ~isunix
        dirs = strread(TempDir,'%s','delimiter',';'); dirs = char(dirs);
    else
        dirs = strread(TempDir,'%s','delimiter',':'); dirs = char(dirs);
    end
    dirs = unique(dirs,'rows');
else
    dirs = InputFolders;
end
Nd = size(dirs,1);
%% Building 4D images
for i = 1:Nd
    disp(['Converting Directory ' num2str(i) ' of ' num2str(Nd)]);
    try
    Niftidir = deblank(dirs(i,:));
    a =dir([Niftidir filesep '*.nii']);
    [namet{1:size(a,1)}] = deal(a.name); TImages = [repmat([Niftidir filesep],[size(a,1) 1]) char(namet)];clear namet;
    if ~isempty(a)&(size(a,1)>1)
        gradient = zeros(size(a,1),3);bvalue = zeros(1,size(a,1));
        [ptht,namt,ext] = fileparts(TImages(1,:));
        [partialnames] = strread(namt,'%s','delimiter','-');
        nname = [char(partialnames{1}) '-' char(partialnames{2}) '-' char(partialnames{3}) '-' char(partialnames{4})];
         Nname = [ptht filesep nname '.nii'];
            for j = 1:size(a,1)
                if exist([Niftidir filesep a(j).name(1:end-4) '.mat'],'file')
                    load([Niftidir filesep a(j).name(1:end-4) '.mat']);
                    gradient(j,:) = ResultsTV.grad;
                    bvalue(j) = ResultsTV.b;
                end
            end
            TImagesdel = TImages;
            indb0 = find(sum(gradient')==0);
            try
                if sum(gradient(:)) ~=0
                    if indb0(end) == size(TImages,1)
                        gradient(indb0(end),:) = [];
                        bvalue(indb0(end)) = [];
                        TImages(indb0(end),:) = [];
                    end
                end
            end
            [Nname] = create_4DVolumes(TImages(1:end,:), Nname);
        if sum(bvalue(:) ~= 0)&sum(gradient(:) ~= 0)
             bvalue = bvalue;
             gradient = gradient';
             save([Nname(1:end-4) '.bvec'],'gradient','-ascii');
             save([Nname(1:end-4) '.bval'],'bvalue','-ascii');
        end
        for j = 1:size(TImagesdel,1)
            delete(deblank(TImagesdel(j,:)));
            if exist([deblank(TImagesdel(j,1:end-4)) '.mat'],'file')
                delete([deblank(TImagesdel(j,1:end-4)) '.mat']);
            end
        end
    elseif ~isempty(a)&(size(a,1)==1)
        [pth,nm,ext] = fileparts(deblank(TImages));
        nmparts = strread(nm,'%s','delimiter','-');
        newname = [char(nmparts{1}) '-' char(nmparts{2}) '-' char(nmparts{3}) '-' char(nmparts{4})];
        movefile(deblank(TImages),[pth filesep newname '.nii']);
        if exist([pth filesep nm '.mat'],'file')
            movefile([pth filesep nm '.mat'],[pth filesep newname '.mat']);
        end
    end
    catch
        deblank(dirs(i,:))
    end
end
return

function [Nname] = create_4DVolumes(Images, Nname);
[pth,name,ext] =fileparts(Images(1,:));
V = spm_vol(Images(1,:));Vt = spm_vol(Images);
for j = 1:size(Images,1)

    Vol(j) = V;
    Vol(j).fname = Nname;
    Vol(j).dt = [spm_type('float32') 0];
    Vol(j).descrip = sprintf('%s','4D Image');
    Vol(j).n = j;
end
Vol =spm_create_vol(Vol);
for j =1:size(Images,1);
    spm_write_vol(Vol(j),spm_read_vols(Vt(j)));
end
Nname = Vol(j).fname;
return
