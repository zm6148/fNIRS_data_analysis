function imgrecon = resetImgRecon(imgrecon)

set(imgrecon.handles.pushbuttonCalcMetrics_new, 'enable','off');
set(imgrecon.handles.editMetricsColormapThreshold_new, 'enable','off');
set(imgrecon.handles.textMetricsColormapThreshold_new, 'enable','off');
if ishandles(imgrecon.handles.hLocalizationError)
    delete(imgrecon.handles.hLocalizationError);
    imgrecon.handles.hLocalizationError = [];
end
if ishandles(imgrecon.handles.hResolution)
    delete(imgrecon.handles.hResolution);
    imgrecon.handles.hResolution = [];
end
imgrecon = setImgReconMetricsColormap(imgrecon, []);
imgrecon.mesh = initMesh();
