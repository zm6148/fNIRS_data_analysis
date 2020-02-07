function probe = pullProbeToHeadsurf(probe,head)

if 1
    probe.optpos_reg = pullPtsToSurf(probe.optpos, head, 'center');
else
    probe.optpos_reg = pullPtsToSurf(probe.optpos, head, 'normal');
end
