function axesv = displayAxesv(axesv, headsurf, digpts)

if isempty(axesv)
    return;
end

zoommore = .8;
if ~digpts.isempty(digpts)
    axesv(1).mode = 'points';
    if headsurf.isempty(headsurf)
        zoommore = 1;
    end
else
    axesv(1).mode = 'surface';
    zoommore = 1.3;
end

view(0,90);
axis(axesv.handles.axesSurfDisplay,'vis3d');
axis(axesv.handles.axesSurfDisplay,'equal');
set(axesv.handles.axesSurfDisplay,'units','normalized');
if strcmp(axesv.mode,'surface')
    pos = get(axesv.handles.axesSurfDisplay, 'position');
    set(axesv.handles.axesSurfDisplay, 'position', [pos(1) pos(2) .5 .4])
end

axesv = setLighting(axesv, headsurf);

set(axesv.handles.axesSurfDisplay, 'cameraupvector', axesv.cameraupvector);

% Usage:      setInitOrientation(axesobj, verr, verl, horr, horl, up, down)
camposvector = [2200.00 2800.00 -2200.00];
switch(headsurf.orientation)
    case 'LIA'
        setOrientation(axesv, 90, 5, 30, 0, 0, 20);
        camposvector = [2200.00 2800.00 -2200.00];
    case 'ASR'
        setOrientation(axesv, 0, 0, 0, 0, 0, 0);
        camposvector = [-2200.00 -2800.00 -2200.00];
    case 'ARI'
        setOrientation(axesv, 0, 0, 0, 0, 0, 0);
        camposvector = [-2200.00 2000.00 2800.00];
end

set(axesv.handles.axesSurfDisplay, {'xlimmode','ylimmode','zlimmode'}, {'manual','manual','manual'});

% Find center of all the objects on the canvas
c1 = [];
c2 = [];
if ~headsurf.isempty(headsurf)
    c1  = headsurf.centerRotation;
end
if ~digpts.isempty(digpts)
    c2  = digpts.center;
end
if ~isempty(c1)
    c = c1;
elseif ~isempty(c2)
    c = c2;
else
    c = [0,0,0];
end

% Find and set camera position 
axesv.campos = c - camposvector;
set(axesv.handles.axesSurfDisplay, 'cameraposition', axesv.campos);

% Set where the camera is pointing to 
if ~all(c==0)
    set(axesv.handles.axesSurfDisplay, 'cameratarget', c);
end

camzoom(axesv.handles.axesSurfDisplay, axesv.zoomincr*zoommore);

