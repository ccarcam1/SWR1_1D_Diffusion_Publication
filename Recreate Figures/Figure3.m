%% Figure 3 (PANEL A)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First open data structure, with variable name "summary"
close all
clearvars -except summary
warning('off') 
figure('Name', 'Panel A')
hold on
counta = 1;
% for i = [3, 7, 12]% this is for the last figure 
% for i = [2, 3, 1]% this is a test of ATP and different salts on SWR1 
for i = [12, 14, 3]
clearvars -except summary i counta slopeval_holding yintercept_holding pvalue_holding rsquared_holding mean_MSD_from_hist x_save y_save
[M, F] = mode(summary(i).line_time);
% Over the course of this project, the line-time changes (this happened
% most frequently at the start, which means I can't average all the traces
% instead i can only average ones with the same line-time, for this reason
% I find the most often occuring line-time in the data set and present the
% mean of this data. 
mode_restriction = summary(i).line_time == M; 
if i ~=12
    fit_restriction = and(summary(i).pvalue<0.05, summary(i).rsquared>0.8);
else
    fit_restriction = mode_restriction;
end
n = sum(and(mode_restriction, fit_restriction));
MSDs_interest = {};
MSDs_interest = summary(i).MSDs(and(mode_restriction, fit_restriction));
mean_MSD_from_hist(counta) = mean(summary(i).dcoef(and(mode_restriction, fit_restriction)));
for loop_var = 1:length(MSDs_interest)
    for j = 1:length(MSDs_interest{loop_var})
        x{j, loop_var} = MSDs_interest{loop_var}(j,1);
        y{j, loop_var} = MSDs_interest{loop_var}(j,2);
    end
end
for loop_var = 1:size(x,1) 
    c = 1;
    for j = 1:size(x,2)
        if ~isempty(x{loop_var,j})
            xx{loop_var, c} = x{loop_var,j};
            yy{loop_var, c} = y{loop_var,j};
            c = c+1;
        end
    end
end
for loop_var = 1:size(xx,1) 
    xxx = [];
    yyy = [];
    c =1;
    for j = 1:size(xx,2)
        if ~isempty(x{loop_var,j})
            xxx(c) = x{loop_var,j};
            yyy(c) = y{loop_var,j};
            c = c+1;
        end
    end
    xm(loop_var) = mean(xxx);
    ym(loop_var) = mean(yyy);
    xs(loop_var) = std(xxx)/sqrt(length(xxx));
    ys(loop_var) = std(yyy)/sqrt(length(yyy));
    xl(loop_var) = length(xxx);
    yl(loop_var) = length(yyy);
end
yl = yl';
    x_save{counta} = xm;
    y_save{counta} = ym;
D = summary(i).meanD;
if i == 3
    coloris = 'g';
elseif i == 7
    coloris = 'b';
elseif i == 12
    coloris = 'k';
else 
    coloris = 'm';
end
    shadedErrorBar(xm, ym, ys, 'lineProps', coloris)
    xm = xm';
    ym = ym';
    fittypeis = strcat('2*', num2str(D), '*x^alpha+c');
    fit_plot = fit(xm(xm<=2), ym(xm<=2), fittype(fittypeis));
    p12 = predint(fit_plot,xm,0.95,'observation','on');
    plot(fit_plot)
    plot(xm,p12,'m--')
xlabel('time')
ylabel('MSD')

xisthis = xm(xm<=2);
yisthis = ym(xm<=2);
semisthis = ys(xm<=2)';

disp(strcat(summary(i).condition, ' has n = ', num2str(n), ' molecules',...
    ' and has an alpha of =', num2str(fit_plot.alpha)))
end
legend off
title('Figure 3: Panel A')
y_upperlim = 0.126258+0.0189;
ylim([0 y_upperlim])
xlim([0 2])

%% (Panel B) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Panel B')

xy = {};
clearvars -except summary
for i = [3, 14, 12]
    if i == 3
        coloris = 'red';
    elseif i == 14
        coloris = 'blue';
    elseif i == 12;
        coloris = 'black';
    end
    tempval = [];
    lengthtime = [];
    trackkept = summary(i).tracking(and(summary(i).bool_kept, summary(i).outlier_exclusion));
    linetimeskept = summary(i).line_time(and(summary(i).bool_kept, summary(i).outlier_exclusion));
    tempval = cellfun(@length, trackkept);
    for j = 1:length(linetimeskept)
    lengthtime(j) = tempval(j)*linetimeskept(j)/1000;
    end
    tracks = {};
    line_time = [];
    tracks = trackkept(lengthtime>8);
    line_time = linetimeskept(lengthtime>8);
    % Seed the random number generator
    rng(4)
    randis = randperm(length(tracks));
    for j = 1:10
    plot(tracks{randis(j)}(:,1)*line_time(randis(j))/1000, (tracks{randis(j)}(:,2)-tracks{randis(j)}(1,2))*0.1, coloris)
    x = []; y = [];
    x = tracks{randis(j)}(:,1)*line_time(randis(j))/1000;
    y = (tracks{randis(j)}(:,2)-tracks{randis(j)}(1,2))*0.1;
    xy{i}{j} = [x(x<8), y(x<8)];
    hold on
%     pause
    end  
    ylim([-2 2])
    xlim([0 8])
end
hold off 
title('Figure 3: Panel B')


%% (Panel C) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(newline)
% Part 1: Exclusion and Outlier Criteria Data refinement
clearvars -except summary
for i = 1:length(summary)
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
    outlier_exclusion = dcoef<numbis;
    summary(i).outlier_exclusion = summary(i).dcoef<numbis;
    summary(i).final_n = length(dcoef(outlier_exclusion));
    summary(i).Median_thresh = median(dcoef(outlier_exclusion));
    summary(i).errorMedian_thresh = (std(dcoef(outlier_exclusion))/sqrt(length(dcoef(outlier_exclusion))))*sqrt(pi/2);

end
forcsv = {};

% each histogram sepately - used Adobe Illustrator to stack 
for i = [14, 3, 13, 4]
    figure('Name', 'Panel C')

    histogram(summary(i).dcoef(and(summary(i).outlier_exclusion, summary(i).bool_kept)),20, 'Normalization', 'probability', 'BinLimits', [0 0.14])
    forcsv{i} = summary(i).dcoef(and(summary(i).outlier_exclusion, summary(i).bool_kept))';
    title('Figure 3: Panel C')
    subtitle(summary(i).condition)
    disp(strcat(summary(i).condition, ' n = ', num2str(sum(and(summary(i).outlier_exclusion, summary(i).bool_kept)))));
end
    
    hold off
%% Panel D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except summary
figure('Name', 'Panel D')

Median_thresh = [summary.Median_thresh].';
errorMedian_thresh = [summary.errorMedian_thresh].';
medianis = Median_thresh;
erroris = errorMedian_thresh;
ATPm = medianis(3);
ATPgsm = medianis(13);
noATPm = medianis(14);
ADPm = medianis(4);
dCas9m =medianis(12);

ATPe = erroris(3);
ATPgse = erroris(13);
noATPe = erroris(14);
ADPe = erroris(4);
dCas9e =erroris(12);


x = [1, 2, 3, 4, 5];
y = [dCas9m,noATPm, ADPm, ATPm,ATPgsm];
SEM = [dCas9e,noATPe, ADPe, ATPe,ATPgse];
errorbar(x, y, SEM, 'ro')
names = {'dCas9';'no ATP';'ADP'; 'ATP'; 'ATP gammaS'; ' '};
set(gca,'xtick',[1:6],'xticklabel',names)
xlim([0 5])
ylim([0 0.05])
title('Panel D')
hold off
%%  Mann-Whitney Median Test 
ATP_i = 3;
ADP_i = 4;
ATPgamms_i = 13;
noATP_i = 14;
dcas9i = 12;
disp(newline)
pairs = [12, 14; 14, 4; 14, 3; 3, 13];
prediction = {'****'; 'n.s.'; '****';'n.s.'};

for i = 1:size(pairs, 1)
one = summary(pairs(i, 1)).dcoef(and(summary(pairs(i, 1)).bool_kept,summary(pairs(i, 1)).outlier_exclusion));%ATP
two = summary(pairs(i, 2)).dcoef(and(summary(pairs(i, 2)).bool_kept,summary(pairs(i, 2)).outlier_exclusion));%ATP
[p,h,stats] = ranksum(one, two);
pis = ['p = ', num2str(p)];
his = ['h = ', num2str(h)];

disp(strcat(summary(pairs(i, 1)).condition, ' and ', summary(pairs(i, 2)).condition, ...
    ' have a null test score of ', his, ' with a significance of: ', pis))
disp(strcat('matching the reported: ',prediction{i}))

end


%% Panel E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Panel E')

% Histcounts to determine which percentage of population is is very slow to stationary phase 
max = 0.14;
min = 0;
bins = 20;
type_of_plot = 'pdf';
type_style = 'stairs';
disp(newline)
for i = 1:length(summary) 
    criteria = and(summary(i).bool_kept, summary(i).bool_kept);
    [summary(i).N,summary(i).edges] = histcounts(summary(i).dcoef(criteria), 10,'BinLimits', [min max]);
    [summary(i).N_immobile, summary(i).edges_immobile] = histcounts(summary(i).dcoef(criteria), bins,'BinLimits', [min max]);
    summary(i).slow = summary(i).N(1);
    summary(i).slow_edge = summary(i).edges(2);
    summary(i).slow100 =  (summary(i).slow/length(summary(i).dcoef(criteria)))*100;
    summary(i).immobile = summary(i).N_immobile(1);
    summary(i).immobile_edge = summary(i).edges_immobile(2);
    summary(i).immobile100 = (summary(i).immobile/length(summary(i).dcoef(criteria)))*100;
disp(strcat(summary(i).condition, ' % immobile is: ', num2str(summary(i).immobile100)))
end
ATP_i = 3;
ADP_i = 4;
ATPgamms_i = 13;
noATP_i = 14;
dcas9i = 12;

% for i = 1:length(summary)
% summary(i).immobile100 = (sum(and(and(summary(i).bool_kept,summary(i).outlier_exclusion), summary(i).dcoef<0.014))/length(summary(i).dcoef(and(summary(i).bool_kept,summary(i).outlier_exclusion))))*100;
% disp(strcat('n = ', num2str(sum(and(summary(i).bool_kept,summary(i).outlier_exclusion)))))
% end
y = [summary(12).immobile100,summary(14).immobile100, summary(4).immobile100,summary(3).immobile100,summary(13).immobile100];

bar(100-y)
names = {'dCas9';'no ATP';'ADP'; 'ATP'; 'ATP gammaS'; ' '};
set(gca,'xtick',[1:6],'xticklabel',names)
ylim([0 100])
title('Panel E')
clearvars -except summary
