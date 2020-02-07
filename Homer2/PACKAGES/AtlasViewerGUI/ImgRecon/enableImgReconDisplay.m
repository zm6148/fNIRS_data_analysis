function enableImgReconDisplay(imgrecon, val)

if val==0
    val='off';
elseif val==1
    val='on';    
end

set(imgrecon.handles.editMetricsColormapThreshold_new, 'enable',val);
set(imgrecon.handles.textMetricsColormapThreshold_new, 'enable',val);
set(imgrecon.handles.editImageReconColormapThreshold, 'enable', val);
set(imgrecon.handles.textImageReconColormapThreshold, 'enable', val);
set(imgrecon.handles.menuItemImageReconGUI, 'enable', val);
