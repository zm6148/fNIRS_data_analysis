function fwmodel = resetSensitivity(fwmodel, probe, dirnameSubj)

% Turn off controls related to projecting fwmodel  onto pialsurf.
% We'll  turn   them back on if probe registration     is successful
set(fwmodel.handles.editSelectChannelSensitivity_new,'enable','off');
set(fwmodel.handles.textSelectChannelSensitivity_new,'enable','off');
set(fwmodel.handles.editSensitivityColormapThreshold_new,'enable','off');
set(fwmodel.handles.textSensitivityColormapThreshold_new,'enable','off');
if ishandles(fwmodel.handles.surf)
    delete(fwmodel.handles.surf);
    fwmodel.handles.surf = [];
end

