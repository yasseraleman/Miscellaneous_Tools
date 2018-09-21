function OutFiles = Select_files(InputDirs,orfilters,andfilters,Ids);
%
% Syntax:
% OutFiles = Select_files(InputDirs,andfilters,Ids);
%
% This function looks(recursively) all files inside the specified folders ("InputDirs"). Only the
% files with the specified search filters in their filenames will be selected. It is
% also possible to select only the files for a specified Subjects Id list.
%
% Input Parameters:
%   InputDir    : Input Search Directories
%   orfilters   : OR Filters. It selects files that contains, at least, one
%                 of the string of the orfilters
%   andfilters  : AND Filters. It selects files that contains all the strings
%                 specified in the variable andfilters.
%   Ids         : Ids lists. It can be a cell array, a char or an Id text
%                 file.
%
% Output Parameters:
%    OutFiles   : Files with the specified characters within the variable "filters" .
%
%
% Example of Usage:
%    OutFiles = Select_files(strvcat('Folder1','Folder2','Folder3'),{'*Connect*','*Dista*'},{'graphmodel'},'IDFilename.txt')
%
%
% See also:
%__________________________________________________________________________
% Authors: Yasser Aleman Gomez
% LIM: HUGGM
% Last update: March 12th 2010
% Version $1.0

warning off
fclose('all');
%% ===================== Checking Input Parameters =======================%
if nargin ==0
    error('Please specify a correct Input directory and correct filters');
    return;
end
if nargin <2
    error('Please specify a correct OR Filter');
    return;
end
%% ===================== End of Checking Input Parameters ================%
%% ============================ Main Program =============================%
Pthst = '';
for i = 1:size(InputDirs,1)
    InputDir = deblank(InputDirs(i,:));
    TempDir = genpath(InputDir);
    if ~isunix
        Pths = strread(TempDir,'%s','delimiter',';'); Pths = char(Pths);
    else
        Pths = strread(TempDir,'%s','delimiter',':'); Pths = char(Pths);
    end
    Pthst = strvcat(Pthst,Pths);
end
pthold = pwd;
% =========================  OR Filtering  ===============================%
OutFiles = '';
cont = 0;
for i = 1:size(Pthst,1)
    %cd(deblank(Pthst(i,:)));
    for j = 1:length(orfilters)
        a = dir([deblank(Pthst(i,:)) filesep char(orfilters{j})]);
        if ~isempty(a)
            [names{1:size(a,1),1}]=deal(a.name);
            files = [repmat([deblank(Pthst(i,:)) filesep],[size(names,1) 1]) char(names)];
            OutFiles = strvcat(OutFiles,files);
            clear names;
        end
    end
end
% =======================  End of OR Filtering   =========================%
% ========================= AND Filtering  ==============================%
if nargin>2
    for k = 1:length(andfilters)
        OutFilest = OutFiles;
        filt = deblank(andfilters{k});
        OutFilest = OutFilest';
        OutFilest =OutFilest(:)';
        ind = strfind(OutFilest,filt);
        [X,Y] = ind2sub(size(OutFiles'),ind);
        rows = unique(Y);
        OutFiles = OutFiles(rows,:);
    end
end
% =======================  End of AND Filtering ==========================%
% ========================= Selecting Ids ================================%
if nargin ==4
    a = whos('Ids');
    switch a.class
        case 'char'
            if exist(Ids,'file')
                Ids = textread(Ids,'%s');
            else
                for i = 1:size(Ids,1)
                    Idst{i} = deblank(Ids(i,:));
                end
                Ids = Idst';
            end
        case 'cell'
            Ids = Ids';
    end
    OutFilest = OutFiles;
    OutFilest = OutFilest';
    OutFilest =OutFilest(:)';
    OutFilesf = '';
    for k = 1:length(Ids)
        filt = deblank(Ids{k});
        ind = strfind(OutFilest,filt);
        if ~isspace(ind)
            [X,Y] = ind2sub(size(OutFiles'),ind);
            rows = unique(Y);
            OutFilesf = strvcat(OutFilesf,OutFiles(rows,:));
        end
    end
    OutFiles = OutFilesf;
end
% ========================= End of Selecting Ids =========================%
%% ========================= End of Main Program =========================%
return