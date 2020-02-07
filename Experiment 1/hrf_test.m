figure(2);
RT = 0.02;
p = [30 35 1 1 6 0 45];
T = 45;
[h,p] = spm_hrf(RT,p,T);
h = h/max(h);
%h = flip(h);
plot(h);

% t=0:1/50:50;
% h2 = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);

hold on;
plot(h2);