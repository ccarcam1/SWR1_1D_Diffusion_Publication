%% Panel A
clearvars -except summary
close all
clc
k = 1;
for i = 1:3
%     grouporder{1, i} = summary(i).condition;
    for j = 1:length(summary(i).dcoef)
        if summary(i).bool_kept(j) == 1
        dcoef(k, 1) = summary(i).dcoef(j); 
        condition{k, 1} = erase(summary(i).condition, ' + 1mM ATP');
        k = k +1;
        end
    end
end

condition = cellstr(condition);
% grouporder = {'25mM KCl + 1mM ATP','70mM KCl + 1mM ATP','200mM KCl + 1mM ATP'};
grouporder = {'25mM KCl','70mM KCl','200mM KCl'};
grouporder = cellstr(grouporder);

figure('Name', 'Panel A')
testmain = violinplot(dcoef(dcoef<0.14), condition(dcoef<0.14),'GroupOrder',grouporder);
testmain(1).ShowMean = 1;
testmain(2).ShowMean = 1;
testmain(3).ShowMean = 1;
title('SWR1')
ylim([0 0.14])

clearvars -except summary
figure('Name', 'Panel A')
k = 1;
for i = 8:10
%     grouporder{1, i} = summary(i).condition;
    for j = 1:length(summary(i).dcoef)
        if summary(i).bool_kept(j) == 1
        dcoef(k, 1) = summary(i).dcoef(j); 
        condition{k, 1} = summary(i).condition;
        k = k +1;
        end
    end
end

condition = cellstr(condition);
grouporder = {'Swc2 25mM noATP','Swc2 70mM noATP','Swc2 150mM noATP'};
grouporder = cellstr(grouporder);
testmain = violinplot(dcoef(dcoef<5), condition(dcoef<5),'GroupOrder',grouporder);
ylim([0 5])
title('Swc2')
testmain(1).ShowMean = 1;
testmain(2).ShowMean = 1;
testmain(3).ShowMean = 1;

%% Panel B
clearvars -except summary
figure('Name', 'Panel B')
len1 = cellfun(@(x) size(x,1), summary(1).tracking).*summary(1).line_time/1000;
len2 = cellfun(@(x) size(x,1), summary(2).tracking).*summary(2).line_time/1000;
len3 = cellfun(@(x) size(x,1), summary(3).tracking).*summary(3).line_time/1000;
Y1 = len1(and(summary(1).bool_kept, summary(1).outlier_exclusion));
Y2 = len2(and(summary(2).bool_kept, summary(2).outlier_exclusion));
Y3 = len3(and(summary(3).bool_kept, summary(3).outlier_exclusion));
[f1, x1] = ecdf(Y1);
[f2, x2] = ecdf(Y2);
[f3, x3] = ecdf(Y3);
typeis = fittype('exp(b*x)');
fit1 = fit(x1-x1(1),1-f1,typeis);
fit2 = fit(x2-x2(1),1-f2,typeis);
fit3 = fit(x3-x3(1),1-f3,typeis);
hold on
s1 = scatter(x1-x1(1),1-f1,5,'filled');
c = get(s1, 'CData');
p1 = plot(fit1);
set(p1, 'Color', c)
s3 = scatter(x3-x3(1),1-f3,5,'filled');
c = get(s3, 'CData');
p3 = plot(fit3);
set(p3, 'Color', c)
s2 = scatter(x2-x2(1),1-f2,5,'filled');
c = get(s2, 'CData');
p2 = plot(fit2);
set(p2, 'Color', c)
t1 = strcat('   t_1_/_2: ',' ', num2str(round((log(2)/-fit1.b),2)), ' s');
t2 = strcat('   t_1_/_2: ',' ', num2str(round((log(2)/-fit2.b),2)), ' s');
t3 = strcat('   t_1_/_2: ',' ', num2str(round((log(2)/-fit3.b),2)), ' s');
legend(strcat(summary(1).condition, ' n = ', num2str(length(Y1))),strcat(' fit', t1), strcat(summary(3).condition, ' n = ', num2str(length(Y3))),strcat(' fit', t3), strcat(summary(2).condition, ' n = ', num2str(length(Y2))),strcat(' fit', t2))
t_half1 = [];
b_interval1 = [];
b_errors1 = [];
t_errors1 = [];
t_half1 = (log(2)/-fit1.b);
b_interval1 = confint(fit1,0.95);
b_errors1 = abs(fit1.b - b_interval1);
t_errors1 = (b_errors1/fit1.b)*t_half1;
t_half2 = (log(2)/-fit2.b);
b_interval2 = confint(fit2,0.95);
b_errors2 = abs(fit2.b - b_interval2);
t_errors2 = (b_errors2/fit2.b)*t_half2;
t_half3 = (log(2)/-fit3.b);
b_interval3 = confint(fit3,0.95);
b_errors3 = abs(fit3.b - b_interval3);
t_errors3 = (b_errors3/fit3.b)*t_half3;
hold off
title('SWR1')
ylabel('1-CDF')
ylim([0 1])
xlabel('dwell-time(s)')
xlim([0 60])

clearvars -except summary
figure('Name', 'Panel B')
len1 = cellfun(@(x) size(x,1), summary(8).tracking).*summary(8).line_time/1000;
len2 = cellfun(@(x) size(x,1), summary(9).tracking).*summary(9).line_time/1000;
len3 = cellfun(@(x) size(x,1), summary(10).tracking).*summary(10).line_time/1000;

Y1 = len1(and(summary(8).bool_kept, summary(8).outlier_exclusion));
Y2 = len2(and(summary(9).bool_kept, summary(9).outlier_exclusion));
Y3 = len3(and(summary(10).bool_kept, summary(10).outlier_exclusion));
[f1, x1] = ecdf(Y1);
[f2, x2] = ecdf(Y2);
[f3, x3] = ecdf(Y3);
typeis = fittype('exp(b*x)');
fit1 = fit(x1-x1(1),1-f1,typeis);
fit2 = fit(x2-x2(1),1-f2,typeis);
fit3 = fit(x3-x3(1),1-f3,typeis);
hold on
s1 = scatter(x1-x1(1),1-f1,5,'filled');
c = get(s1, 'CData');
p1 = plot(fit1);
set(p1, 'Color', c)
s3 = scatter(x3-x3(1),1-f3,5,'filled');
c = get(s3, 'CData');
p3 = plot(fit3);
set(p3, 'Color', c)
s2 = scatter(x2-x2(1),1-f2,5,'filled');
c = get(s2, 'CData');
p2 = plot(fit2);
set(p2, 'Color', c)
t1 = strcat('   t_1_/_2: ',' ', num2str(round((log(2)/-fit1.b),2)), ' s');
t2 = strcat('   t_1_/_2: ',' ', num2str(round((log(2)/-fit2.b),2)), ' s');
t3 = strcat('   t_1_/_2: ',' ', num2str(round((log(2)/-fit3.b),2)), ' s');
legend(strcat(summary(1).condition, ' n = ', num2str(length(Y1))),strcat(' fit', t1), strcat(summary(3).condition, ' n = ', num2str(length(Y3))),strcat(' fit', t3), strcat(summary(2).condition, ' n = ', num2str(length(Y2))),strcat(' fit', t2))
t_half1 = [];
b_interval1 = [];
b_errors1 = [];
t_errors1 = [];
t_half1 = (log(2)/-fit1.b);
b_interval1 = confint(fit1,0.95);
b_errors1 = abs(fit1.b - b_interval1);
t_errors1 = (b_errors1/fit1.b)*t_half1;
t_half2 = (log(2)/-fit2.b);
b_interval2 = confint(fit2,0.95);
b_errors2 = abs(fit2.b - b_interval2);
t_errors2 = (b_errors2/fit2.b)*t_half2;
t_half3 = (log(2)/-fit3.b);
b_interval3 = confint(fit3,0.95);
b_errors3 = abs(fit3.b - b_interval3);
t_errors3 = (b_errors3/fit3.b)*t_half3;
hold off
title('Swc2')
ylabel('1-CDF')
ylim([0 1])
xlabel('dwell-time(s)')
xlim([0 15])

%% Panel C
clearvars -except summary
kk = 1;
for i = [1,2,3,4,12,13,14,7,8,9,10]
    if any([8,9,10] ==i)
        numbis = 5;
    else
        numbis = 0.14;
    end
    dcoef = [];
    k = 1;
    for j = 1:length(summary(i).dcoef)
        if summary(i).bool_kept(j) == 1
        dcoef(k) = summary(i).dcoef(j);
        k = k+1;
        end
    end
    SWR1_stats(kk).condition = summary(i).condition; 
    SWR1_stats(kk).decoef = dcoef;
    SWR1_stats(kk).oldnval = length(dcoef); 
    SWR1_stats(kk).newnval = length(dcoef(dcoef<numbis));
    SWR1_stats(kk).Median_thresh = median(dcoef(dcoef<numbis));
    SWR1_stats(kk).errorMedian_thresh = (std(dcoef(dcoef<numbis))/sqrt(length(dcoef(dcoef<numbis))))*sqrt(pi/2);
    SWR1_stats(kk).meanD_oldexcl = mean(dcoef); 
    SWR1_stats(kk).stdmean_oldexcl = std(dcoef); 
    SWR1_stats(kk).semmean_oldexcl = std(dcoef)/sqrt(length(dcoef)); 
    kk = kk+1;
end

figure('Name', 'Panel C')
mdc_ATP_SWR1200 = SWR1_stats(1).Median_thresh;
SEMdc_ATP_SWR1200 = SWR1_stats(1).errorMedian_thresh;
dcoef_ATP_SWR1200 = SWR1_stats(1).decoef(SWR1_stats(1).decoef<0.14);
whatisit = 2;
mdc_ATP_SWR125 = SWR1_stats(whatisit).Median_thresh;
SEMdc_ATP_SWR125 = SWR1_stats(whatisit).errorMedian_thresh;
dcoef_ATP_SWR125 = SWR1_stats(whatisit).decoef(SWR1_stats(whatisit).decoef<0.14);
whatisit = 3;
mdc_ATP_SWR1 = SWR1_stats(whatisit).Median_thresh;
SEMdc_ATP_SWR1 = SWR1_stats(whatisit).errorMedian_thresh;
dcoef_ATP_SWR1 = SWR1_stats(whatisit).decoef(SWR1_stats(whatisit).decoef<0.14);
whatisit = 7;
mdc_noATP_SWR1 = SWR1_stats(whatisit).Median_thresh;
SEMdc_noATP_SWR1 = SWR1_stats(whatisit).errorMedian_thresh;
dcoef_noATP_SWR1 = SWR1_stats(whatisit).decoef(SWR1_stats(whatisit).decoef<0.14);
whatisit = 4;
mdc_ADP_SWR1 = SWR1_stats(whatisit).Median_thresh;
SEMdc_ADP_SWR1 = SWR1_stats(whatisit).errorMedian_thresh;
dcoef_ADP_SWR1 = SWR1_stats(whatisit).decoef(SWR1_stats(whatisit).decoef<0.14);
whatisit = 6;
mdc_ATPgS_SWR1 = SWR1_stats(whatisit).Median_thresh;
SEMdc_ATPgS_SWR1 = SWR1_stats(whatisit).errorMedian_thresh;
dcoef_ATPgS_SWR1 = SWR1_stats(whatisit).decoef(SWR1_stats(whatisit).decoef<0.14);
whatisit = 5;
mdc_dCas9 = SWR1_stats(whatisit).Median_thresh;
SEMdc_dCas9 = SWR1_stats(whatisit).errorMedian_thresh;
dcoef_dCas9 = SWR1_stats(whatisit).decoef(SWR1_stats(whatisit).decoef<0.14);
whatisit = 9;
mdc_Swc2150 = SWR1_stats(whatisit).Median_thresh;
SEMdc_Swc2150 = SWR1_stats(whatisit).errorMedian_thresh;
dcoef_Swc2150 = SWR1_stats(whatisit).decoef(SWR1_stats(whatisit).decoef<5);
whatisit = 10;
mdc_Swc225 = SWR1_stats(whatisit).Median_thresh; 
SEMdc_Swc225 = SWR1_stats(whatisit).errorMedian_thresh;
dcoef_Swc225 = SWR1_stats(whatisit).decoef(SWR1_stats(whatisit).decoef<5);
whatisit = 11;
mdc_Swc270 = SWR1_stats(whatisit).Median_thresh;
SEMdc_Swc270 =SWR1_stats(whatisit).errorMedian_thresh;
dcoef_Swc270 = SWR1_stats(whatisit).decoef(SWR1_stats(whatisit).decoef<5);

x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
y = [mdc_dCas9,mdc_noATP_SWR1, mdc_ADP_SWR1,mdc_ATP_SWR125,mdc_ATP_SWR1, mdc_ATP_SWR1200,mdc_ATPgS_SWR1,mdc_Swc225,mdc_Swc270, mdc_Swc2150];
SEM = [SEMdc_dCas9,SEMdc_noATP_SWR1, SEMdc_ADP_SWR1,SEMdc_ATP_SWR125,SEMdc_ATP_SWR1,SEMdc_ATP_SWR1200,SEMdc_ATPgS_SWR1,SEMdc_Swc225,SEMdc_Swc270,SEMdc_Swc2150];
hold on
x_dcoef = ones(1, length(dcoef_dCas9));
scatter(x_dcoef.*1, dcoef_dCas9, '_');
clear x_dcoef 
x_dcoef = ones(1, length(dcoef_ATP_SWR1200));
scatter(x_dcoef.*6, dcoef_ATP_SWR1200, '_');
clear x_dcoef  
x_dcoef = ones(1, length(dcoef_ATP_SWR125));
scatter(x_dcoef.*4, dcoef_ATP_SWR125, '_');
clear x_dcoef 
x_dcoef = ones(1,length(dcoef_ATP_SWR1));
scatter(x_dcoef.*5, dcoef_ATP_SWR1, '_');
clear x_dcoef 
x_dcoef = ones(1, length(dcoef_noATP_SWR1));
scatter(x_dcoef.*2, dcoef_noATP_SWR1, '_');
clear x_dcoef 
x_dcoef = ones(1, length(dcoef_ADP_SWR1));
scatter(x_dcoef.*3, dcoef_ADP_SWR1, '_');
clear x_dcoef 
x_dcoef = ones(1, length(dcoef_ATPgS_SWR1));
scatter(x_dcoef.*7, dcoef_ATPgS_SWR1, '_');
clear x_dcoef 
x_dcoef = ones(1, length(dcoef_Swc225));
scatter(x_dcoef.*8, dcoef_Swc225, '_');
clear x_dcoef 
x_dcoef = ones(1, length(dcoef_Swc270));
scatter(x_dcoef.*9, dcoef_Swc270, '_');
clear x_dcoef 
x_dcoef = ones(1, length(dcoef_Swc2150));
scatter(x_dcoef.*10, dcoef_Swc2150, '_');
clear x_dcoef
errorbar(x, y, SEM, 'ko')
names = {'dCas9';'no ATP';'ADP'; 'SWR1 ATP 25mM KCl';'SWR1 ATP 70mM KCl'; 'SWR1 ATP 200mM KCl' ; 'ATP gammaS' ; 'Swc2 25mM';'Swc2 70mM'; 'Swc2 150mM'};
set(gca,'xtick',[1:10],'xticklabel',names)
uncoupled_lim_SWR1 = 3.67E+01; %um2/sec
coupled_lim_SWR1_small = 1.05E-01; % um2/sec
coupled_lim_SWR1_large = 1.83E-01; %um2/sec
yline([uncoupled_lim_SWR1 coupled_lim_SWR1_small coupled_lim_SWR1_large],'--',{'uncoupled','coupled(min)', 'coupled(max)'})
set(gca,'yscale','log')
uncoupled_lim_Swc2 = 1.25E+02; %um2/sec
coupled_lim_Swc2_small = 4.01E+00; % um2/sec
coupled_lim_Swc2_large = 6.86E+00; %um2/sec
yline([uncoupled_lim_Swc2 coupled_lim_Swc2_small coupled_lim_Swc2_large],'--',{'uncoupled','coupled(min)', 'coupled(max)'})
set(gca,'yscale','log')
title('Theoretical Upper Limits to Diffusion SWR1 and Swc2')