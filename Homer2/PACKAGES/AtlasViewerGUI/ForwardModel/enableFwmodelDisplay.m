function fwmodel = enableFwmodelDisplay(fwmodel, val)

if val==0
    val='off';
elseif val==1
    val=val;
end

set(fwmodel.handles.popupmenuImageDisplay,'enable',val);
set(fwmodel.handles.editSelectChannelSensitivity_new,'enable',val);
set(fwmodel.handles.textSelectChannelSensitivity_new,'enable',val);
set(fwmodel.handles.editSensitivityColormapThreshold_new,'enable',val);
set(fwmodel.handles.textSensitivityColormapThreshold_new,'enable',val);
set(fwmodel.handles.menuItemGetSensitivityatMNICoordinates,'enable',val);
