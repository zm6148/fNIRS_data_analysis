function status = headvol2headsurf(dirname, headvol)

status=1;

if ~exist('dirname','var')
    dirname='./';
elseif dirname(end)~='/' && dirname(end)~='\'
    dirname(end+1)='/';
end

if ~exist('headvol','var')
    headvol=[];
end

headvol = getHeadvol(headvol, dirname);
if isempty(headvol.img)
    return;
end

fv = isosurface(headvol.img,.9);

% isosurface flips x and y, so we have to either flip x and y back, or have
% the transform do it.
fv.vertices=[fv.vertices(:,2) fv.vertices(:,1) fv.vertices(:,3)];

% Since we flipped mesh x and y, no need for transformation 
headsurf2vol = eye(4);

h = waitbar(0,'Downsampling head surface. This may take a few minutes...');
[fv.vertices,fv.faces] = meshresample(fv.vertices,fv.faces,0.15);
close(h);

if ~exist([dirname 'anatomical'], 'dir')
    mkdir([dirname 'anatomical']);
end
write_surf([dirname 'anatomical/headsurf.mesh'], fv.vertices, fv.faces);
save([dirname 'anatomical/headsurf2vol.txt'],'-ascii','headsurf2vol');

status=0;