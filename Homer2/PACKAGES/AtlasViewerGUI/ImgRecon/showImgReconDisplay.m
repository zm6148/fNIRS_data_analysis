function imgrecon = showImgReconDisplay(imgrecon, hAxes, valLocErr, valRes, valHbO, valHbR)

if isempty(hAxes)
    hAxes=gca;
end

if strcmp(valRes, 'on') | strcmp(valLocErr, 'on')
    valMetrics='on';
else
    valMetrics='off';
end
if strcmp(valHbO, 'on') | strcmp(valHbR, 'on')
    valImgRecon='on';
else
    valImgRecon='off';
end


if strcmp(valMetrics,'off')
    imgrecon = setImgReconMetricsColormap(imgrecon, []);
end
if strcmp(valImgRecon,'off')
    imgrecon = setImgReconColormap(imgrecon, []);
end

if strcmp(valMetrics,'on')
    imgrecon = setImgReconMetricsColormap(imgrecon, hAxes);    
end
if strcmp(valImgRecon,'on')
    imgrecon = setImgReconColormap(imgrecon, hAxes);    
end

set(imgrecon.handles.hLocalizationError, 'visible',valLocErr);
set(imgrecon.handles.hResolution, 'visible',valRes);
set(imgrecon.handles.hHbO, 'visible',valHbO);
set(imgrecon.handles.hHbR, 'visible',valHbR);

set(imgrecon.handles.editImageReconColormapThreshold, 'visible',valImgRecon);
set(imgrecon.handles.textImageReconColormapThreshold, 'visible',valImgRecon);
set(imgrecon.handles.editMetricsColormapThreshold_new, 'visible',valMetrics);
set(imgrecon.handles.textMetricsColormapThreshold_new, 'visible',valMetrics);


