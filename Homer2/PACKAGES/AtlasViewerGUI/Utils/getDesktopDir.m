function [desktopDir, err] = getDesktopDir()

desktopDir = '';
err=true;
if ispc()
    [status, username] = dos('echo %USERNAME%');
    if ~isalnum(username(end))
        username(end)='';
    end
    desktopDir = sprintf('c:/users/%s/desktop/', username);
    if ~exist(desktopDir, 'dir')~=7
         desktopDir = 'c:/users/public/';
    end
elseif isunix() | ismac()
    if exist('~/Desktop/','dir')
        currdir = pwd;
        if ~exist('~/Desktop','file')
            mkdir('~/Desktop/');
        end
        cd('~/Desktop/');
        desktopDir = pwd;
        cd(currdir);
    end
end
if exist(desktopDir, 'dir')~=7
    try
        mkdir(desktopDir);
        err = false;
    catch
        err = true;
    end
else
    err = false;
end
if err
    desktopDir='';
elseif desktopDir(end) ~= '/'
    desktopDir(end+1) = '/';
end




% ---------------------------------------------------------
function y = isalnum( x )

y = all( (x>='0' & x<='9') | (x>='a' & x<='z') | (x>='A' & x<='Z') );


