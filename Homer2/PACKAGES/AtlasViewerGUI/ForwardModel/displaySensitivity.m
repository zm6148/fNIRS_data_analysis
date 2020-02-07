function fwmodel = displaySensitivity(fwmodel, pialsurf, labelssurf, probe, hAxes)

if isempty(fwmodel)
    return;
end
if fwmodel.isempty(fwmodel)
    set(fwmodel.handles.menuItemImageReconGUI,'enable','off');
    return;
end

val = get(fwmodel.handles.popupmenuImageDisplay,'value');

if ~exist('hAxes','var')
    hAxes = fwmodel.handles.axes;
end

if leftRightFlipped(fwmodel)
    axes_order = [2,1,3];
else
    axes_order = [1,2,3];
end

if ~isempty(probe.optpos_reg)
    set(fwmodel.handles.menuItemGenerateMCInput,'enable','on');
else
    set(fwmodel.handles.menuItemGenerateMCInput,'enable','off');    
end

% Error checks
if isempty(probe.optpos_reg)
    return;
end
if isempty(probe.ml)
    return;
end
iCh = find(probe.ml(:,1)==fwmodel.Ch(1) & probe.ml(:,2)==fwmodel.Ch(2), 1);
if fwmodel.Ch(1)==0 & fwmodel.Ch(2)==0
    iCh = 0;
end
if isempty(iCh)
    return;
end
set(fwmodel.handles.menuItemImageReconGUI,'enable','on');

if ishandles(fwmodel.handles.surf)
    delete(fwmodel.handles.surf);
end

% Wavelength to display always one for now. TBD: Add feature to select
% between wavelengths
iW = 1;

viewOrigin(hAxes);
if iCh>0
    if all(fwmodel.Adot(iCh(1),:,iW)==0)
        intensity = fwmodel.cmThreshold(1).*ones(size(fwmodel.Adot(iCh(1),:,iW),2), 1);
    else
        intensity = log10(fwmodel.Adot(iCh(1),:,iW));
    end
else
    if all(fwmodel.Adot(:,:,iW)==0)
        intensity = fwmodel.cmThreshold(1).*ones(size(sum(fwmodel.Adot(:,:,iW),1),2), 1);
    else
        intensity = log10(sum(fwmodel.Adot(:,:,iW), 1));
    end
end

fwmodel.handles.surf = ....
    displayIntensityOnMesh(fwmodel.mesh, intensity, fwmodel.cmThreshold, ...
                           fwmodel.colormin, 'off','off', axes_order);
hold off;


if ishandles(fwmodel.handles.surf)
    
    fwmodel = showFwmodelDisplay(fwmodel, hAxes, 'on');
    fwmodel = enableFwmodelDisplay(fwmodel, 'on');
    if val~=1
        set(fwmodel.handles.popupmenuImageDisplay,'value',1);
    end
    fwmodel = setSensitivityColormap(fwmodel, hAxes);
    
    probe = setProbeDisplay(probe);
    set(probe.handles.hMeasList,'color','y');
    if iCh>0
        set(probe.handles.hMeasList(iCh(1),1),'color','g','linewidth',2);
        set(probe.handles.hMeasList(iCh(1),1),'visible','on');
    else
        set(probe.handles.hMeasList(:,1),'color','g','linewidth',2);
        set(probe.handles.hMeasList(:,1),'visible','on');
    end
    
    set(pialsurf.handles.radiobuttonShowPial, 'value',0);
    uipanelBrainDisplay(pialsurf.handles.radiobuttonShowPial, {pialsurf, labelssurf});
else
    
    fwmodel = enableFwmodelDisplay(fwmodel, 'off');    
    fwmodel = setSensitivityColormap(fwmodel, []);
    
end

