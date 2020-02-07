function probe = resetProbeGui(probe)

hSprings = probe.handles.hSprings;
hOptodes = probe.handles.hOptodes;
hOptodesCircles = probe.handles.hOptodesCircles;
hMeasList = probe.handles.hMeasList;
hMeasCortex = probe.handles.hMeasCortex;
hRays = probe.handles.hRays;

if ishandles(hOptodes)
   delete(hOptodes);
end
if ishandles(hOptodesCircles)
   delete(hOptodesCircles);
end
if ishandles(hMeasCortex)
   delete(hMeasCortex);
end
if ishandles(hMeasList)
   delete(hMeasList);
end
if ishandles(hSprings)
   delete(hSprings);
end  
if ishandles(hRays)
   delete(hRays);
end  


