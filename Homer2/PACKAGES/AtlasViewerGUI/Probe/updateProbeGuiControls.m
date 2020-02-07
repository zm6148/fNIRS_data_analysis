function probe = updateProbeGuiControls(probe,headsurf,method)

if ~exist('method','var') || isempty(method)
    if ~ishandles(hSprings)
        method = 'digpts';
    else
        method = 'springs';
    end
end

if strcmp(method,'springs') 
    set(probe.handles.checkboxHideProbe,'enable','on');
    set(probe.handles.checkboxHideSprings,'enable','on');
    set(probe.handles.checkboxHideDummyOpts,'enable','on');
    set(probe.handles.checkboxHideMeasList,'enable','on');
elseif strcmp(method,'digpts')
    set(probe.handles.checkboxHideProbe,'enable','on');
    set(probe.handles.checkboxHideSprings,'enable','off');
    set(probe.handles.checkboxHideDummyOpts,'enable','off');
    set(probe.handles.checkboxHideMeasList,'enable','on');
end

% Figure out if the probe is pre-registered. That it is either in 
% position to be pulled toward the head or if we can register it using 
% springs and anchor points, if they exist. The button 
% pushbuttonRegisterProbeToSurface is used for both cases but has to be 
% enabled
if exist('headsurf','var')
    if ~isempty(headsurf)
        p = probe.optpos;
        [~, ~, d] = nearest_point(headsurf.mesh.vertices, p);
        
        % Check if proble is flat and has 
        if ~isempty(probe.sl) & ~isempty(probe.al)
            b = true;
        elseif all(d<30) & std(d,1,2)<10
            b = true;
        else 
            b = false;
        end
    else
        b = false;
    end
else
    b = true;
end

if ~isempty(probe.optpos) & b==1
    set(probe.handles.pushbuttonRegisterProbeToSurface,'enable','on');
else
    set(probe.handles.pushbuttonRegisterProbeToSurface,'enable','off');
end
if ~isempty(probe.optpos)
    set(probe.handles.checkboxOptodeSDMode,'enable','on');
    set(probe.handles.checkboxOptodeCircles,'enable','on');
else
    set(probe.handles.checkboxOptodeSDMode,'enable','off');
    set(probe.handles.checkboxOptodeCircles,'enable','off');
end
if ~isempty(probe.ml)
    set(probe.handles.checkboxHideMeasList,'enable','on');
else
    set(probe.handles.checkboxHideMeasList,'enable','off');
end
if ~isempty(probe.sl)
    set(probe.handles.checkboxHideSprings,'enable','on');
    set(probe.handles.checkboxHideDummyOpts,'enable','on');
else
    set(probe.handles.checkboxHideSprings,'enable','off');
    set(probe.handles.checkboxHideDummyOpts,'enable','off');
end
if ~isempty(probe.optpos_reg)
    set(probe.handles.menuItemSaveRegisteredProbe,'enable','on');
    set(probe.handles.menuItemProjectOptodesToCortex, 'enable','on');
    set(probe.handles.menuItemProjectChannelsToCortex, 'enable','on');
else
    set(probe.handles.menuItemSaveRegisteredProbe,'enable','off');
    set(probe.handles.menuItemProjectOptodesToCortex, 'enable','off');
    set(probe.handles.menuItemProjectChannelsToCortex, 'enable','off');
end

probe.hideProbe     = get(probe.handles.checkboxHideProbe,'value');
probe.hideSprings   = get(probe.handles.checkboxHideSprings,'value');
probe.hideDummyOpts = get(probe.handles.checkboxHideDummyOpts,'value');
probe.hideMeasList  = get(probe.handles.checkboxHideMeasList,'value');
