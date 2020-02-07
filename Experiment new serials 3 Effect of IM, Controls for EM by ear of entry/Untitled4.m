
d_new = [];
for ii = 1:40
    d_new = cat(2, d_new, -d(:,ii));
end

d = d_new;

save('TOT_3_20191101_1032_03.mat');