function fwmodel = showFwmodelDisplay(fwmodel, hAxes, val)

if isempty(hAxes)
    hAxes=gca;
end

if val==0
    val='off';
elseif val==1
    val=val;
end

if strcmp(val,'off')
    fwmodel = setSensitivityColormap(fwmodel, []);
elseif strcmp(val,'on')
    fwmodel = setSensitivityColormap(fwmodel, hAxes);
end

set(fwmodel.handles.editSelectChannelSensitivity_new,'visible',val);
set(fwmodel.handles.textSelectChannelSensitivity_new,'visible',val);
set(fwmodel.handles.editSensitivityColormapThreshold_new,'visible',val);
set(fwmodel.handles.textSensitivityColormapThreshold_new,'visible',val);
set(fwmodel.handles.surf,'visible',val);


