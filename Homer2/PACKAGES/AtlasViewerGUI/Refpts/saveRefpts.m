function refpts = saveRefpts(refpts, dirname, T_vol2mc, mode)

if isempty(refpts) | isempty(refpts.pos)
    return;
end

if ~exist('dirname','var')
    dirname = [refpts.pathname, '/anatomical/'];
else
    if dirname(end)~='/' && dirname(end)~='\'
        dirname(end+1)='/';
end
    dirname = [dirname, '/anatomical/'];    
end
if ~exist(dirname, 'dir')
    mkdir(dirname);
end
if ~exist('T_vol2mc','var') | isempty(T_vol2mc)
    T_vol2mc = eye(4);
end
if ~exist('mode', 'var')
    mode='nooverwrite';
end

if ~isempty(refpts.handles.selected)
    for ii=1:size(refpts.pos,1)
        if ishandles(refpts.handles.selected(ii))
            delete(refpts.handles.selected(ii));
            refpts.handles.selected(ii) = -1;
        end
    end
end

if strcmp(mode,'nosave')
    return;
end

if exist([dirname,'refpts.txt'],'file') & strcmp(mode,'overwrite')
    menu(sprintf('Old refpts.txt were moved to refpts.txt.bak'),'OK');
    movefile([dirname, 'refpts.txt'], [dirname, 'refpts.txt.bak']);
    if exist([dirname,'refpts_labels.txt'],'file')
        movefile([dirname, 'refpts_labels.txt'], [dirname, 'refpts_labels.txt.bak']);
    end
    if exist([dirname,'refpts2vol.txt'],'file')
        movefile([dirname, 'refpts2vol.txt'], [dirname, 'refpts2vol.txt.bak']);
    end
end

if ~exist([dirname 'refpts2vol.txt'], 'file') | strcmp(mode, 'overwrite')
    T_2vol = refpts.T_2vol;
    save([dirname 'refpts2vol.txt'], 'T_2vol', '-ascii');
end

if ~exist([dirname 'refpts.txt'], 'file') | strcmp(mode, 'overwrite')
    % Unapply T_2vol to get back to original ref pts
    pos = refpts.pos;
    pos = xform_apply(pos, inv(T_vol2mc * T_2vol));
    save([dirname 'refpts.txt'], 'pos', '-ascii');
end

if ~exist([dirname 'refpts_labels.txt'], 'file') | strcmp(mode, 'overwrite')
    fd = fopen([dirname 'refpts_labels.txt'], 'wt');
    for ii=1:length(refpts.labels)
        fprintf(fd, '%s\n', refpts.labels{ii});
    end
    fclose(fd);
end


