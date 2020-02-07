function probe = getProbe(probe, dirname, headsurf, displaymode)

if iscell(dirname)
    for ii=1:length(dirname)
        probe = getProbe(probe, dirname{ii}, headsurf, displaymode);
        if ~probe.isempty(probe)
            return;
        end
    end
    return;
end

if isempty(dirname)
    return;
end

if exist([dirname, 'probe.mat'],'file')

    load([dirname, 'probe.mat'],'-mat');
    
else
  
    filetype = 0;
    inputfile='';
    if exist(dirname,'file')==2
        inputfile = dirname;
        k = find(dirname=='/' | dirname=='\');
        dirname = dirname(1:k(end));
        filetype=2;
    elseif exist(dirname,'dir')==7
        if dirname(end)~='/' && dirname(end)~='\'
            dirname(end+1)='/';
        end
        if exist([dirname 'probe.txt'],'file')
            inputfile = [dirname, 'probe.txt'];
            filetype=1;
        elseif exist([dirname 'digpts.txt'],'file')
            inputfile = [dirname, 'digpts.txt'];
            filetype=2;
        else
            return;
        end 
    else
        return;
    end
    
    if filetype==2
        digpts = getDigpts([],inputfile);
        probe.optpos   = [digpts.srcpos; digpts.detpos];
        probe.nsrc     = size(digpts.srcpos,1);
        probe.noptorig = size(probe.optpos,1);
        probe.ndet     = probe.noptorig - probe.nsrc;
        probe.srcmap   = digpts.srcmap;
        probe.detmap   = digpts.detmap;
        probe.orientation = digpts.orientation; 
        probe.center   = digpts.center; 
        
        % Get measurement list. However the measurement list is only 
        % valid if srcmap and detmap are 1:nsrc and 1:ndet respectively.
        probe = loadMeasdList(probe, dirname);
        
    elseif filetype==1
    
        probe.optpos = load(inputfile,'-ascii');
        probe.srcpos = probe.optpos;
        probe.nsrc = size(probe.srcpos, 1);
        probe.ndet = 0;
        
    end
end
    



% -----------------------------------------
function b = isnumber(str)    
b = ~isempty(str2num(str));




% -----------------------------------------
function msg = warningMsgs(filenm, probe, SD)

msg = '';
n1 = size(SD.SrcPos, 1);
n2 = size(SD.DetPos, 1);
m1 = probe.nsrc;
m2 = probe.ndet;

msg1 = sprintf('Warning: The number of optodes in digpts.txt (srcs: %d, dets: %d)\n', m1, m2); 
msg2 = sprintf('does not equal the number in the SD file %s (srcs: %d, dets: %d). ', filenm, n1, n2);

msg = [msg1, msg2];



% -----------------------------------------
function probe = loadMeasdList(probe, dirname)

if probe.noptorig<=0
    return;
end
if ~all((1:probe.nsrc)==probe.srcmap)
    return;
end
if ~all((1:probe.ndet)==probe.detmap)
    return;
end

foos1 = dir([dirname '*.nirs']);
foos2 = dir([dirname '*.SD']);
foos3 = dir([dirname '*.sd']);

files = [foos1; foos2; foos3];

for ii=1:length(files)
    s = load(files(ii).name,'-mat');
    if size(s.SD.SrcPos,1)==probe.nsrc & size(s.SD.DetPos,1)==probe.ndet
        SD = s.SD;
        break;
    end
end
if ~exist('SD','var')
    return;
end

smax = max(SD.MeasList(:,1));
dmax = max(SD.MeasList(:,2));
if smax>probe.nsrc | dmax>probe.ndet
    msg1 = sprintf('Warning: Measurement list indices in the SD file %s exceed\n', files(1).name);
    msg2 = sprintf('the number of optodes. Will load only valid channels.');
    menu([msg1, msg2], 'OK');
    drawnow;    
    ks = find(SD.MeasList(:,1) > probe.nsrc);
    kd = find(SD.MeasList(:,2) > probe.ndet);
    SD.MeasList([ks,kd],:) = [];
end

if ~isempty(SD.MeasList)
    k = SD.MeasList(:,4)==1;
    probe.ml = SD.MeasList(k,:);
end

if ~probe.isempty(probe)
    probe.pathname = dirname;
end
