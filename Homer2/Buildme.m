function Buildme(appName, dependList, exclList)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Buildme allows building of the current directories project in 
% without having to change the output directory or update the .m 
% files list, every time (as seems to be the case with deploytool)
%
% Also it finds all the .m files under the current directory 
%

DEBUG = 1;

if isunix() || ismac()
    disp(sprintf('Buildme is not yet implmented for this platform.'));
end

currDir = pwd;

% Error check args
if ~exist('appName','var') 
    [pp,fs] = getpathparts(currDir);
    appName = pp{end};
end
if ~exist('dependList','var')
    dependList = {};
end
if ~exist('exclList','var')
    exclList = {};
end

% Get all .m files which will go into making the app executable
appFiles = findDotMFiles('.',exclList);
for ii=1:length(dependList)
    appFiles1 = findDotMFiles(dependList{ii},exclList);
    appFiles = [appFiles, appFiles1];
end

appFileMain = sprintf('%s.m', appName);
targetName = sprintf('.%s%s.exe', filesep, appName);


%%% Go through all the apps, contruct a string listing all the .m files 
%%% on which each app depends and then compile the app using mcc.
compileSwitches = '-w enable:specified_file_mismatch -w enable:repeated_file -w enable:switch_ignored -w enable:missing_lib_sentinel -w enable:demo_license';
appFilesStr = sprintf('-v %s%s%s', currDir, filesep, appFileMain');

% Remove main m file and remove Buildme.m from app files list 
appFiles = removeEntryFromList(appFileMain, appFiles);
appFiles = removeEntryFromList('Buildme.m', appFiles);

% Construct files list portion of build command
if DEBUG
    fid = fopen('Buildme.log','w');
end
for jj=2:length(appFiles)
    appFilesStr = sprintf('%s -a %s', appFilesStr, appFiles{jj});
    if DEBUG
        fprintf(fid, '%s\n', appFiles{jj});
    end
end

% Complete the final build command and execute it
buildcmdstr = sprintf('mcc -o %s -W main:%s -T link:exe -d %s %s %s', appName, appName, currDir, compileSwitches, appFilesStr);
disp(buildcmdstr);
eval(buildcmdstr);




% -------------------------------------------------------------------------
function dotmfiles = findDotMFiles(subdir, exclList)

if ~exist('exclList','var')
    exclList = {};
end

dotmfiles = {};
currdir = pwd;

if exist(subdir, 'dir')~=7
    fprintf('Warning: folder %s doesn''t exist under %s\n', subdir, pwd);
    return;
end
cd(subdir);

% If current subjdir is in the exclList then go back to curr dir and exit
subdirFullpath = pwd;

for ii=1:length(exclList)
    if ~isempty(findstr(exclList{ii}, subdirFullpath))
        cd(currdir);
        return;
    end
end

files = dir('*');
if isempty(files)
    return;
end

for ii=1:length(files)
    exclFlag = false;
    if isdotmfile(files(ii))
        for kk=1:length(exclList)
            if strcmp(files(ii).name, exclList{kk})
                exclFlag = true;
            end
        end
        if exclFlag==true
            continue;
        end
        dotmfiles{end+1} = sprintf('%s%s%s', pwd, filesep, files(ii).name);
    elseif files(ii).isdir && ~iscurrdir(files(ii)) && ~isparentdir(files(ii))
        dotmfiles = [dotmfiles, findDotMFiles(files(ii).name, exclList)];
    end
end
cd(currdir);



% -------------------------------------------------------------------------
function b = isdotmfile(file)

b=0;
if file.isdir
    return;
end
if file.name(end) ~= 'm' || file.name(end-1) ~= '.'
    return;
end
b=1;



% -------------------------------------------------------------------------
function b = iscurrdir(file)

b=0;
if ~file.isdir
    return;
end
if isempty(file.name)
    return;
end
if isempty(file.name)
    return;
end
if length(file.name)==1
    if file.name(1)~='.'
        return;
    end
end
if (length(file.name)==2)
    if (file.name(1)~='.') || (file.name(2)~='/' && file.name(2)~='\')
        return;
    end
end
if (length(file.name)>2)
    return;
end

b=1;



% -------------------------------------------------------------------------
function b = isparentdir(file)

b=0;
if ~file.isdir
    return;
end
if isempty(file.name)
    return;
end
if isempty(file.name)
    return;
end
if length(file.name)==1
    return;
end
if (length(file.name)==2)
    if (file.name(1)~='.') || (file.name(2)~='.')
        return;
    end
end
if (length(file.name)==3)
    if (file.name(1)~='.') || (file.name(2)~='.') || (file.name(2)~='/' && file.name(2)~='\')
        return;
    end
end
if (length(file.name)>3)
    return;
end
b=1;



% -------------------------------------------------------------------------
% Helper function: remove name arg from list
function list = removeEntryFromList(name, list)

temp = strfind(list, name);
k=[];
for ii=1:length(temp)
    if ~isempty(temp{ii})
        k=ii;
    end
end
list(k) = [];
