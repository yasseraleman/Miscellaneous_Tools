function OutFiles = sel_files(InputDir,filtro);
TempDir = genpath(InputDir);

if ~isunix
    Pthst = strread(TempDir,'%s','delimiter',';'); Pthst = char(Pthst);
else
    Pthst = strread(TempDir,'%s','delimiter',':'); Pthst = char(Pthst);
end
pthold = pwd;
OutFiles = '';
cont = 0;
for i = 1:size(Pthst,1)
    cd(deblank(Pthst(i,:)));
    a = dir(filtro);
    if ~isempty(a)
    [names{1:size(a,1),1}]=deal(a.name);
    files = [repmat([deblank(Pthst(i,:)) filesep],[size(names,1) 1]) char(names)];
    OutFiles = strvcat(OutFiles,files);
    end
end
cd(pthold);
return