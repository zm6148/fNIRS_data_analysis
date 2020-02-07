function ml=getMeasListFromAxes(optpos_src, optpos_det, h_edges)

ml=zeros(length(h_edges),2);
for i=1:length(h_edges)
    x=get(h_edges(i),'xdata');
    y=get(h_edges(i),'ydata');
    z=get(h_edges(i),'zdata');
    p1=[x(1) y(1) z(1)];
    p2=[x(2) y(2) z(2)];
    [foo j]=nearest_point(optpos_src,p1);
    ml(i,1)=j;
    [foo k]=nearest_point(optpos_det,p2);
    ml(i,2)=k;
end
