function imgrecon = setImgReconColormap(imgrecon, hAxes)

hclim = imgrecon.handles.editImageReconColormapThreshold;

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
if ~ishandles(imgrecon.handles.colorbar)
    imgrecon.handles.colorbar = colorbar;
    set(imgrecon.handles.colorbar, 'visible','off', 'position',[0.68, 0.45, 0.02, 0.30]);
else
    imgrecon.cmThreshold = cmThreshold;
end
    
val = get(imgrecon.handles.popupmenuImageDisplay,'value');
if (val==4 | val==5) & ~isempty(hAxes)
    set(imgrecon.handles.colorbar, 'visible','on');
    cm = jet(100);
    colormap(cm);
    caxis(cmThreshold);
else
    set(imgrecon.handles.colorbar, 'visible','off');
end

