function fwmodel = setSensitivityColormap(fwmodel, hAxes)

hclim = fwmodel.handles.editSensitivityColormapThreshold_new;

if ~isempty(hAxes)
    axes(hAxes);
end

% Error checking
clim = caxis;
cmThreshold = str2num(get(hclim,'string'));
if isempty(cmThreshold)
    set(hclim,'string',sprintf('%g %g',clim(1),clim(2)));
    return;
end
if length(cmThreshold)~=2
    set(hclim,'string',sprintf('%g %g',clim(1),clim(2)));
    return;
end
if cmThreshold(1)>=cmThreshold(2)
    set(hclim,'string',sprintf('%g %g',clim(1),clim(2)));
    return;
end

% Passed error checks. Now set the colorbar colormap
if ~ishandles(fwmodel.handles.colorbar)
    fwmodel.handles.colorbar = colorbar;
    set(fwmodel.handles.colorbar, 'visible','off', 'position',[0.68, 0.45, 0.02, 0.30]);
else
    fwmodel.cmThreshold = cmThreshold;
end

% Visible or not based on what image popupmenu option was selected
val = get(fwmodel.handles.popupmenuImageDisplay,'value');
if val==1 & ~isempty(hAxes)
    set(fwmodel.handles.colorbar, 'visible','on');
    cm = jet(100);
    cm(1,:) = fwmodel.colormin;
    colormap(cm);
    caxis(cmThreshold);    
else
    set(fwmodel.handles.colorbar, 'visible','off');
end


