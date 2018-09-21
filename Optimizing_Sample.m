function [OutIdsFile, OutDemoFile ] = Optimizing_Sample(IdsFile, DemoFile);
%
% Syntax :
%  [OutIdsFile, OutDemoFile ] = Optimizing_Sample(IdsFile, DemoFile);
%
% This scripts extracts a paired sample from an unpaired sample. 
%
% Input Parameters:
%       IdsFile                 : Unpaired Ids
%       DemoFile                : Demographic File for unpaired subjects
%
% Output Parameters:
%       OutIdsFile              : Paired Ids
%       OutDemoFile             : Paired Demographic File
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% December 1st 2014
% Version $1.0


IdsFile = '/media/Data/PEH/IDs.txt';
DemoFile = '/media/Data/PEH/hsp_1.txt';

% Creating output filenames
[pthid, nmid, extid] = fileparts(IdsFile);
OutIdsFile = [pthid filesep nmid '_paired.txt'];

[pthde, nmde, extde] = fileparts(DemoFile);
OutDemoFile = [pthde filesep nmde '_paired.txt'];

% Reading Ids File
[Ids]=textread(IdsFile,'%s');
SubjIds = char(Ids);

% Reading Demographic Data
[group,Sex,Age] = textread(DemoFile,'%u %u %f','delimiter',' ','headerlines',1);

indp = find(group == 1); % Patients must be idenfitied by 1
Idsp = SubjIds(indp,:);
indc = find(group == 0); % Controls must be idenfitied by 0
Idsc = SubjIds(indc,:);

Sexp = Sex(indp);
Sexc = Sex(indc);
Agep = Age(indp);
Agec = Age(indc);
indfp = find(Sexp == 1); % Females must be idenfitied by 1
indmp = find(Sexp == 2); % Males must be idenfitied by 2
indfc = find(Sexc == 1);
indmc = find(Sexc == 2);

% Females
Np = length(indfp);
Nc = length(indfc);
totalDistm = 10^12;
totalDistf = 10^12;
Nperm = 1000;
for k = 1:Nperm
    disp(['Permutation Number: ' num2str(k)]);
    clear indtemp indt;
    % Females
    Np = length(indfp);
    Nc = length(indfc);
    indrp = randperm(Np);
    indrc = randperm(Nc);
    samplep = Agep(indfp(indrp),:);
    samplec = Agec(indfc(indrc),:);
    idsPatients = Idsp(indfp(indrp),:);
    idsControls = Idsc(indfc(indrc),:);
    if Np>=Nc
        sample = pdist2(samplep,samplec);
    else
        sample = pdist2(samplec,samplep);
    end
    [M,N] = size(sample);
    sumDist = 0;
    for i = 1:N
        S = sample(:,i);
        
        if i == 1
            ind = find(S == min(S));
            indt(i) = ind(1);
        else
            ind = find(ismember([1:M],indt)==0);
            indtemp = find(S(ind) == min(S(ind)));
            indt(i) = ind(indtemp(1));
        end
        sumDist = sumDist + sample(indt(i),i);
    end
    if totalDistf>=sumDist
        totalDistf = sumDist;
        if Np>=Nc
            sampleControlsF = samplec;
            samplePatientsF = samplep(indt,:);
            idsControlFinalF = idsControls;
            idsPatientsFinalF = idsPatients(indt,:);
        else
            sampleControlsF = samplec(indt,:);
            samplePatientsF = samplep;
            idsControlFinalF = idsControls(indt,:);
            idsPatientsFinalF = idsPatients;
        end
    end
    
    % Males
    clear indtemp indt;
    Np = length(indmp);
    Nc = length(indmc);
    indrp = randperm(Np);
    indrc = randperm(Nc);
    samplep = Agep(indmp(indrp),:);
    samplec = Agec(indmc(indrc),:);
    idsPatients = Idsp(indmp(indrp),:);
    idsControls = Idsc(indmc(indrc),:);
    if Np>=Nc
        sample = pdist2(samplep,samplec);
    else
        sample = pdist2(samplec,samplep);
    end
    [M,N] = size(sample);
    sumDist = 0;
    for i = 1:N
        S = sample(:,i);
        
        if i == 1
            ind = find(S == min(S));
            indt(i) = ind(1);
        else
            ind = find(ismember([1:M],indt)==0);
            indtemp = find(S(ind) == min(S(ind)));
            indt(i) = ind(indtemp(1));
        end
        sumDist = sumDist + sample(indt(i),i);
    end
    if totalDistm>=sumDist
        totalDistm = sumDist;
        if Np>=Nc
            sampleControlsM = samplec;
            samplePatientsM = samplep(indt,:);
            idsControlFinalM = idsControls;
            idsPatientsFinalM = idsPatients(indt,:);
        else
            sampleControlsM = samplec(indt,:);
            samplePatientsM = samplep;
            idsControlFinalM = idsControls(indt,:);
            idsPatientsFinalM = idsPatients;
        end
    end
end
IdsFinale = strvcat(idsControlFinalF,idsControlFinalM,idsPatientsFinalF,idsPatientsFinalM);
groups = [[sampleControlsF(:,1);sampleControlsM(:,1)]*0;[samplePatientsF(:,1);samplePatientsM(:,1)]*0+1];
Sex = [sampleControlsF(:,1)*0+1;sampleControlsM(:,1)*0+2;samplePatientsF(:,1)*0+1;samplePatientsM(:,1)*0+2];
Ages = [sampleControlsF;sampleControlsM;samplePatientsF;samplePatientsM];

cad = [repmat('     ',[size(Sex,1) 1]) num2str(groups) repmat('     ',[size(Sex,1) 1]) num2str(Sex) repmat('     ',[size(Sex,1) 1]) num2str(Ages)];

% Saving Ids for paired sample
fid = fopen(OutIdsFile,'wt');
N = size(IdsFinale,1);
for i = 1:N
    fprintf(fid,'%s\n',deblank(IdsFinale(i,:)));
end
fclose(fid);

% Saving Demographic Data for paired sample
fid = fopen(DemoFile,'wt');
N = size(cad,1);
for i = 1:N
    fprintf(fid,'%s\n',deblank(cad(i,:)));
end
fclose(fid);
return;

