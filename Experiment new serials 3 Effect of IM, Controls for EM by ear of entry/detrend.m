function dod_derend = detrend( dod, order )
% 20th poly detrend fucntion

dod_derend = [];

for ii = 1 : size(dod,2)
    
    data = dod(:,ii);
    t = (1:length(data))';
    
    [p,s,mu] = polyfit(t,data,order);
    f_y = polyval(p,t,[],mu);
    
    dt_data = data - f_y;
    
    dod_derend = cat(2, dod_derend, dt_data);

end

