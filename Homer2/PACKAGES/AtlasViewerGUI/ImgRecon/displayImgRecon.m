function imgrecon = displayImgRecon(imgrecon, fwmodel, pialsurf, labelssurf, probe, hAxes)

if isempty(imgrecon)
    return;
end
if fwmodel.isempty(fwmodel)
    return;
end

% Since sensitivity profile exists, enable all image panel controls
% for calculating metrics
set(imgrecon.handles.pushbuttonCalcMetrics_new, 'enable','on');

if imgrecon.isempty(imgrecon)
    return;
end

if leftRightFlipped(imgrecon)
    axes_order = [2,1,3];
else
    axes_order = [1,2,3];
end


val = get(imgrecon.handles.popupmenuImageDisplay,'value');

if ~exist('hAxes','var')
    hAxes = imgrecon.handles.axes;
end

% Error checks
if isempty(imgrecon.localizationError)
    return;
end
if isempty(imgrecon.resolution)
    return;
end
if isempty(probe.optpos_reg)
    return;
end
if isempty(probe.ml)
    return;
end

viewOrigin(hAxes);
intensityLC  = imgrecon.localizationError;
intensityRes = imgrecon.resolution;
if ~isempty(imgrecon.mesh)
    imgrecon.handles.hLocalizationError = ....
        displayIntensityOnMesh(imgrecon.mesh, intensityLC, imgrecon.cmThresholdMetrics, ...
                               imgrecon.colormin, 'off','off', axes_order);
    imgrecon.handles.hResolution = ....
        displayIntensityOnMesh(imgrecon.mesh, intensityRes, imgrecon.cmThresholdMetrics, ...
                               imgrecon.colormin, 'off','off', axes_order);
end
hold off;


% Enable of disable display controls based on the availability of the 
% hLocalizationError or hResolution handles
if ishandles(imgrecon.handles.hLocalizationError) 
    
    set(imgrecon.handles.editMetricsColormapThreshold_new, 'enable','on');
    set(imgrecon.handles.textMetricsColormapThreshold_new, 'enable','on');
    
    % Turn sensitivity display off
    fwmodel = showFwmodelDisplay(fwmodel, hAxes, 'off');
    
    % Turn localization error on and resolution display off
    imgrecon = showImgReconDisplay(imgrecon, hAxes, 'on', 'off', 'off', 'off');
    
    if val~=2
        set(imgrecon.handles.popupmenuImageDisplay,'value',2);
    end
    
    set(pialsurf.handles.radiobuttonShowPial, 'value',0);
    uipanelBrainDisplay(pialsurf.handles.radiobuttonShowPial, {pialsurf, labelssurf});
elseif val==2

    set(imgrecon.handles.editMetricsColormapThreshold_new, 'enable','off');
    set(imgrecon.handles.textMetricsColormapThreshold_new, 'enable','off');
   
    % Turn localization error on and resolution display off
    imgrecon = showImgReconDisplay(imgrecon, hAxes, 'off', 'off', 'off', 'off');
    
end


% Enable of disable display controls based on the availability of the 
% hLocalizationError or hResolution handles
if ishandles(imgrecon.handles.hLocalizationError)
    set(imgrecon.handles.editMetricsColormapThreshold_new, 'enable','on');
    set(imgrecon.handles.textMetricsColormapThreshold_new, 'enable','on');

    % Turn sensitivity display off
    fwmodel = showFwmodelDisplay(fwmodel, hAxes, 'off');
    
    % Turn localization error on and resolution display off
    imgrecon = showImgReconDisplay(imgrecon, hAxes, 'off', 'on', 'off', 'off');    
    
    if val~=3
        set(imgrecon.handles.popupmenuImageDisplay,'value',3);
    end
    
    set(pialsurf.handles.radiobuttonShowPial, 'value',0);
    uipanelBrainDisplay(pialsurf.handles.radiobuttonShowPial, {pialsurf, labelssurf});
else

    set(imgrecon.handles.editMetricsColormapThreshold_new, 'enable','off');
    set(imgrecon.handles.textMetricsColormapThreshold_new, 'enable','off');
   
    % Turn localization error on and resolution display off
    imgrecon = showImgReconDisplay(imgrecon, hAxes, 'off', 'off', 'off', 'off');
      
end

