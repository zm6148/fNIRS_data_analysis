function threshold = set_threshold(optpos)

dm = distmatrix(optpos);
dmin = min(dm(find(dm(:)~=0)));
if isempty(dmin)
    threshold(1)=1;
    threshold(2)=1;
    threshold(3)=1;    
else
    threshold(1)=dmin/3.75;
    threshold(2)=dmin/3.40;
    threshold(3)=dmin/3.25;
end