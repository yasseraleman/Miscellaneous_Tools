function Outfiles = Extract_brain(T1Images,T2Images);
%
% Syntax :
% Outfiles = Extract_brain(T1Images,T2Images);
%
% This function computes  Brain masks for MRI data. It creates also the GM,
% WM and CSF segmentation files previously obtained by a segmentation process.
%
% Input Parameters:
%  T1Images       : T1 Weighted Image.
%  T2Images       : T2 Weighted Image.
%
%
% Output Parameter:
%   Outfiles       : Brain mask, GM, WM and CSF segmentation files
%
% See also: i_spm_segment
%__________________________________________________
% Authors: Yasser Aleman Gomez
%LIM
% February 1th 2012
% Version $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempd = which('spm_preproc8.m');
Outfiles = '';
if ~isempty(tempd)
    V = spm_vol(T1Images);Vs = spm_vol(T2Images);
    Ns = size(V,1);
    for i = 1:Ns
        T1 = deblank(T1Images(i,:));
        %outnames = change_spacen(Vs(i).fname,V(i).fname,7);Vt = spm_vol(outnames);
        [path,name,ext] = fileparts(T1);
        job.data{1} = T1;
        job.channel(1).vols{1} = V.fname;
        job.channel(1).biasreg = 1.0000e-04;
        job.channel(1).biasfwhm = 60;
        job.channel(1).write = [0 0];
        
%         job.channel(2).vols{1} = Vs.fname;
%         job.channel(2).biasreg = 1.0000e-04;
%         job.channel(2).biasfwhm = 60;
%         job.channel(2).write = [0 0];
        
        job.opts.tpm = {[which('TPM.nii')]};
        job.opts.ngaus = [2 2 2 3 4 2];
        job.opts.biasreg = 1.0000e-04;
        job.opts.biasfwhm = 60;
        job.opts.affreg = 'mni';
        job.opts.warpreg = 4;
        job.opts.samp = 3;
        job.extopts.dartelwarp = 1;
        job.extopts.sanlm = 1;
        job.extopts.mrf = 0.1500;
        job.extopts.cleanup = 1;
        job.extopts.print = 1;
        job.output.GM.native = 1;
        job.output.GM.warped = 0;
        job.output.GM.modulated = 0;
        job.output.GM.dartel = 0;
        job.output.WM.native = 1;
        job.output.WM.warped = 0;
        job.output.WM.modulated = 0;
        job.output.WM.dartel = 0;
        job.output.CSF.native = 1;
        job.output.CSF.warped = 0;
        job.output.CSF.modulated = 0;
        job.output.CSF.dartel = 0;
        job.output.bias.native = 0;
        job.output.bias.warped = 0;
        job.output.bias.affine = 0;
        job.output.label.native = 1;
        job.output.label.warped = 0;
        job.output.label.dartel = 0;
        job.output.jacobian.warped = 0;
        job.output.warps = [0 0];
        job.bias = [0 0 0];
        job.label = [1 0 0 0];
        job.jacobian = 0;
        job.biasreg = 1.0000e-04;
        job.biasfwhm = 60;
        job.warp.affreg = 'mni';
        job.warp.samp = 3;
        job.warp.reg = 4;
        job.warp.write = [0 0];
        job.warp.sanlm = 1;
        job.warp.mrf = 0.1500;
        job.warp.print = 1;
        job.warp.cleanup = 1;
        job.warp.dartelwarp = 1;
        job.warps = [0 0];
        native = [1 0 0;1 0 0;1 0 0;0 0 0;0 0 0;0 0 0];
        for j = 1:6
            job.tissue(j).tpm = {[which('TPM.nii') ',' num2str(i)]};
            job.tissue(j).ngaus = job.opts.ngaus(j);
            job.tissue(j).native = native(j,:);
            job.tissue(j).warped = [0 0 0];
        end
        cg_vbm8_run(job);
        Outfiles = strvcat(Outfiles,[repmat([path filesep],[4 1]) strvcat('p0','p1','p2','p3') repmat([name ext(1:4)],[4 1])]);
        Outfiles = strvcat(Outfiles,[path filesep name '_seg8.mat'],[path filesep 'p' name '_seg8.txt']);
%         v1 = spm_vol([path filesep 'T1.nii']);
%         v2 = spm_vol([path filesep 'p0T1.nii']);
%         v1_1 = spm_read_vols(v1);
%         v2_1 = spm_read_vols(v2);
%         v2_1 = logical(v2_1).*v1_1;
%         spm_write_vol(v1, v2_1);
    end
else
    global defaults
    spm_defaults;
    dseg				= defaults.segment;
    V = spm_vol(T1Images);Vs = spm_vol(T2Images);
    Ns = size(V,1);
    VG0 = which('T1.nii');
    for i =1:Ns
        disp(['Extracting    ' num2str(i) ' off ' num2str(Ns) ': ======>>   ' deblank(T1Images(i,:))]);
        outnames = change_spacen(Vs(i).fname,V(i).fname,7);Vt = spm_vol(outnames);
        i_spm_segment(V(i),VG0,dseg,Vt);
        [pth,nm,ext] = fileparts(V(i).fname);
        [pth1,nm1,ext1] = fileparts(Vt.fname);
        GMImages = [pth filesep 'p1' nm ext(1:4)];
        WMImages = [pth filesep 'p2' nm ext(1:4)];
        CSFImages = [pth filesep 'p3' nm ext(1:4)];
        Create_Brain_Masks(GMImages, WMImages, CSFImages);
        delete(outnames);
        %delete([pth filesep 'p1' nm ext(1:4)]);
        %delete([pth filesep 'p2' nm ext(1:4)]);
        %delete([pth filesep 'p3' nm ext(1:4)]);
        delete([pth filesep 'p4' nm ext(1:4)]);
        delete([pth filesep 'p5' nm ext(1:4)]);
        delete([pth filesep 'p6' nm ext(1:4)]);
        [pth1, nm1,ext1] = fileparts(GMImages);
        Outfiles = strvcat(Outfiles,[pth1 filesep nm1(3:end) '_BM' ext1(1:4)],GMImages,WMImages,CSFImages,[pth filesep nm '_seg8.mat'],[pth1 filesep 'y_' nm ext(1:4)]);
        %delete([pth filesep nm '_seg8.mat']);
        %delete([pth1 filesep 'y_' nm1 ext1(1:4)]);
    end
end
return