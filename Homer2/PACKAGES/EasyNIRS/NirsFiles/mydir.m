function dirs = mydir(str)

if ~exist('str','var')
    dirs = dir();
else
    dirs = dir(str);
end

if isempty(dirs)
    fields = fieldnames(dirs);
    args = '';
    for ii=1:length(fields)
        if ii==1
            args = [args '''' fields{ii} ''',' '{}'];
        else
            args = [args ',''' fields{ii} ''',' '{}'];
        end
    end
    args = [ args ',''idx'',0,''subjdir'',{},''subjdiridx'',0,''filename'',{}' ];
    eval(sprintf('dirs = struct(%s);',args));
else
    for jj=1:length(dirs)
        dirs(jj).idx = 0;
        dirs(jj).subjdir = '';                    
        dirs(jj).subjdiridx = 0;
        dirs(jj).filename = dirs(jj).name;
    end 
end