function probe = initProbeProjection(probe)

% dynamic handles
if ishandles(probe.handles.hMeasCortex)
   delete(probe.handles.hMeasCortex);
   probe.handles.hMeasCortex=[];
end
if ishandles(probe.handles.hMeasToLabelsProjTbl)
   delete(probe.handles.hMeasToLabelsProjTbl);
   probe.handles.hMeasToLabelsProjTbl=[];
end
if ishandles(probe.handles.hRays)
   delete(probe.handles.hRays);
   probe.handles.hRays=[];
end


% static handles
if isempty(probe.optpos_reg)
    set(probe.handles.menuItemProjectOptodesToCortex, 'enable','off');
    set(probe.handles.menuItemProjectChannelsToCortex, 'enable','off');    
else
    set(probe.handles.menuItemProjectOptodesToCortex, 'enable','on');
    set(probe.handles.menuItemProjectChannelsToCortex, 'enable','on');    
end

probe.ptsProj_cortex = [];
probe.ptsProj_cortex_mni = [];

