function [nz,iz,rpa,lpa,cz] = getLandmarks(refpts)

nz  = [];
iz  = [];
rpa = [];
lpa = [];
cz  = [];

if isempty(refpts)
    return;
end

for ii=1:length(refpts.labels)
    switch(lower(refpts.labels{ii}))
        case 'nz'
            nz  = refpts.pos(ii,:);
        case 'iz'
            iz  = refpts.pos(ii,:);
        case {'rpa','a2','ar'}
            rpa = refpts.pos(ii,:);
        case {'lpa','a1','al'}
            lpa = refpts.pos(ii,:);
        case 'cz'
            cz  = refpts.pos(ii,:);
    end
end
