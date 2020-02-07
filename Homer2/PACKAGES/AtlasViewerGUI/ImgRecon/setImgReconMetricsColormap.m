function imgrecon = setImgReconMetricsColormap(imgrecon, hAxes)

hclim = imgrecon.handles.editMetricsColormapThreshold_new;

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

% Visible or not based on what iamge popup option was selected
val = get(imgrecon.handles.popupmenuImageDisplay,'value');
if (val==2 | val==3) & ~isempty(hAxes)
    set(imgrecon.handles.colorbar, 'visible','on');
    cm = jet(100);
    cm(1,:) = imgrecon.colormin;
    colormap(cm);
    caxis(cmThreshold);    
else
    set(imgrecon.handles.colorbar, 'visible','off');
end

