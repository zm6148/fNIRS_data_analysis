function fluenceProf = loadFluenceProf(filenm, varargin)

fluenceProf = [];

% Load fluence profile from file
fields = varargin;
if isempty(fields)
    fprintf('Load command: load(''%s'')\n', filenm);
    fluenceProf0 = load(filenm);
else
    strLst=[];
    n=length(fields);
    for ii=1:n
        if ii<n
            strLst = [strLst '''' fields{ii} ''','];
        elseif ii==n
            strLst = [strLst '''' fields{ii} ''''];
        end
    end    
    cmd = sprintf('fluenceProf0 = load(''%s'', %s, ''-mat'');', filenm, strLst);
    % fprintf('Load command: %s\n', cmd);
    eval(cmd);
end


% Make sure the loaded profile follows current format
fluenceProf = initFluenceProf();
fields = fieldnames(fluenceProf);
n=length(fields);
for ii=1:n
    if isfield(fluenceProf0, fields{ii})
        eval(sprintf('fluenceProf.%s = fluenceProf0.%s;',fields{ii},fields{ii}));
    end
end


