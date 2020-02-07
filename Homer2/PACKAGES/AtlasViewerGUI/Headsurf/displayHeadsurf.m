function headsurf = displayHeadsurf(headsurf, hAxes)

if isempty(headsurf)
    return;
end
if headsurf.isempty(headsurf)
    return;
end
if ~exist('hAxes','var')
    hAxes = headsurf.handles.axes;
end

if ishandles(headsurf.handles.surf)
    delete(headsurf.handles.surf);
end

if leftRightFlipped(headsurf)
    axes_order = [2,1,3];
else
    axes_order = [1,2,3];
end

if isempty(headsurf.center)
    headsurf.center = findcenter(fv.vertices);
    headsurf.centerRotation = headsurf.center;
end
if isempty(headsurf.centerRotation)
    headsurf.centerRotation = headsurf.center;
end

h=[];
if isempty(headsurf.mesh)
    menu('head file does not exist in current directory','ok');
    return;
else
    viewOrigin(hAxes);
    h = viewsurf(headsurf.mesh, .7, headsurf.color, 'off', axes_order);
    hold off
end
headsurf.handles.surf = h;

if ishandles(headsurf.handles.surf)
    if ishandles(headsurf.handles.radiobuttonShowHead)
        set(headsurf.handles.radiobuttonShowHead,'value',1);
        set(headsurf.handles.radiobuttonShowHead,'enable','on');
        set(headsurf.handles.editTransparency,'enable','on');
        set(headsurf.handles.editTransparency,'string', num2str(get(headsurf.handles.surf,'facealpha')));
        set(headsurf.handles.menuItemMakeProbe,'enable','on');
        set(headsurf.handles.menuItemImportProbe,'enable','on');
    end
else
    if ishandles(headsurf.handles.radiobuttonShowHead)
        set(headsurf.handles.radiobuttonShowHead,'enable','off');
        set(headsurf.handles.menuItemMakeProbe,'enable','off');
        set(headsurf.handles.menuItemImportProbe,'enable','off');
    end
end


