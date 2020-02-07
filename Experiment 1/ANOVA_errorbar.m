% ranova error bar
 a = [0.2484271867	0.2373349511
0.01011778774	-0.1950746709
0.9305528312	0.706446293
0.1661238585	-0.1052209826
-0.5878643357	-0.6983923424
0.02160390771	-0.07666936779
0.03832340679	0.1759634648
1.2745707	1.346419589
0.2200663947	0.3997211342
-0.09502402451	-0.8290842423
0.02108435057	-0.0944659964
-0.07032529152	0.03940802387
-0.2575138622	-0.1408736929
-0.7290626996	-0.8804283855];



b= [0.2744546011	0.3751194472
0.2693644001	0.1860099726
0.6912706046	0.5537235553
0.05591398668	-0.03562716219
0.09437134301	-0.05466049785
0.3366712231	0.4880253002
0.2156257957	0.3016170072
2.126252108	1.42372314
0.02909973038	0.3231645312
1.030929285	0.1111516677
0.4723284584	0.2589618362
-0.02722171773	-0.337733469
0.2365958446	-0.01902711567
-0.04529256165	-0.2252122426];
% figure;bar(mean(b))

data_by_sub = a';
subj = size(data_by_sub,2);
condition_mean = mean(data_by_sub,2);
condition_std_error  = std(data_by_sub,0,2)./sqrt(subj);

diff_from_mean = kron(ones(size(data_by_sub,1),1),mean(data_by_sub,1)) - kron(ones(size(data_by_sub)),mean(mean(data_by_sub,1),2));

adjusted_data = data_by_sub - diff_from_mean;
adjusted_std_error = std(adjusted_data,0,2)./sqrt(subj);

figure;bar([1,2],mean(a))
hold on;errorbar([1,2],mean(a),adjusted_std_error,'linestyle','none')


data_by_sub = b';
subj = size(data_by_sub,2);
condition_mean = mean(data_by_sub,2);
condition_std_error  = std(data_by_sub,0,2)./sqrt(subj);

diff_from_mean = kron(ones(size(data_by_sub,1),1),mean(data_by_sub,1)) - kron(ones(size(data_by_sub)),mean(mean(data_by_sub,1),2));

adjusted_data = data_by_sub - diff_from_mean;
adjusted_std_error = std(adjusted_data,0,2)./sqrt(subj);

hold on;bar([4,5],mean(b))
hold on;errorbar([4,5],mean(b),adjusted_std_error,'linestyle','none')


set(gca,'xtick',[1.5,4.5],'xticklabel',{'right tgPCS';'right cIFS'})