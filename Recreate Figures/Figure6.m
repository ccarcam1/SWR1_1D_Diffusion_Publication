%% Figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find compaction ratio 
clc
clearvars -except data summary
close all
k = 1;
k2 = 1;
k3 = 1;
unwrappings = [];
fivepN_unwrappings = [];
fivepN_all = [];

for i = 1:length(data)
    if not(matches(data(i).condition, 'No Nucleosome Control'))
        for j = 1:length(data(i).type)
            if and(j == 1, ~isempty(data(i).unwrappingcounts2))
                unwrappings(k) = data(i).unwrappingcounts2;
                fivepN_unwrappings(k)= data(i).fivepNdis(1);
                fivepN_all(i) = data(i).fivepNdis(1);
                i_is(k) = i;
                k = k+1;
            elseif j == 1
                fivepN_all(i) = data(i).fivepNdis(1);
            elseif j == 2
                fivepN_nuc_ret(k3) = data(i).fivepNdis(2); 
                k3 = k3+1;
            end
        end
    else
        for j = 1:length(data(i).type)
            if j == 1
                fivepN_nonucs_ext(k2) = data(i).fivepNdis(1);
            elseif j == 2
                fivepN_nonucs_ret(k2) = data(i).fivepNdis(2);
                k2 = k2+1;
            end
        end
    end
end


%% determine the total number of nucleosomes/ array
% close all
figure('Name', 'Supplementary Figure')
scatter(fivepN_unwrappings(fivepN_unwrappings>12), unwrappings(fivepN_unwrappings>12))
[P, S] = polyfit(fivepN_unwrappings(fivepN_unwrappings>12), unwrappings(fivepN_unwrappings>12),1);
polyval(P, fivepN_unwrappings(fivepN_unwrappings>12));
yfit = P(1)*fivepN_unwrappings(fivepN_unwrappings>12)+P(2);

alpha = 0.05; % Significance level
[yfit,delta] = polyconf(P,fivepN_unwrappings(fivepN_unwrappings>12),S,'alpha',alpha);


hold on;
plot(fivepN_unwrappings(fivepN_unwrappings>12),yfit-delta,'r--',fivepN_unwrappings(fivepN_unwrappings>12),yfit+delta,'r--')
plot(fivepN_unwrappings(fivepN_unwrappings>12),yfit,'r-.');
distance_per_nuc = -1/P(1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yresid = unwrappings(fivepN_unwrappings>12) - yfit;
% Square the residuals and total them to obtain the residual sum of squares:
SSresid = sum(yresid.^2);
% Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
SStotal = (length(unwrappings(fivepN_unwrappings>12))-1) * var(unwrappings(fivepN_unwrappings>12));
% Compute R2 using the formula given in the introduction of this topic:
rsq = 1 - SSresid/SStotal; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('Supplementary Figure')
disp(strcat('the distance per nucleosome is: ', num2str(distance_per_nuc)))
disp(strcat('the r2 for the fit is: ', num2str(rsq)))
disp(strcat('n = ', num2str(sum(fivepN_unwrappings>12))))
disp(strcat('mx+b:', num2str(P(1)), 'x +', num2str(P(2))))

%%
figure('Name', 'Main Figure')
% histogram(unwrappings(fivepN_unwrappings>12), 8)
% hold on
histogram(round((mean(fivepN_nuc_ret)-fivepN_all)/distance_per_nuc), 8)
hold off

mean_nucleosome_count = mean(round((mean(fivepN_nuc_ret)-fivepN_all)/distance_per_nuc));
SEM_nucleosome_count = std(round((mean(fivepN_nuc_ret)-fivepN_all)/distance_per_nuc))/sqrt(length(round((mean(fivepN_nuc_ret)-fivepN_all)/distance_per_nuc)));
mean_lengthDNA = mean(fivepN_unwrappings);

mean_distance_bw_nucs = mean_lengthDNA/mean_nucleosome_count;
disp(newline)
disp(strcat('The mean number of nucs is: ', num2str(mean_nucleosome_count)))
disp(strcat('The SEM is : ', num2str(SEM_nucleosome_count)))

disp(strcat('The mean distance between ucleosomes is: ', num2str(mean_distance_bw_nucs)))
disp(strcat('n = ', num2str(length(round((mean(fivepN_nuc_ret)-fivepN_all)/distance_per_nuc)))))
copypaste = [];
copypaste = round((mean(fivepN_nuc_ret)-fivepN_all)/distance_per_nuc)';

%% (Panel A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Panel A')
hold on
n_is_unwrap = 0;
n_is_retrac = 0;
n_is_refer = 0;
for i = 1:length(data)
%     if i == 10
%         line_width = 3;
%     else
        line_width = 0.5;
%     end
    if matches(data(i).condition, 'Nucleosome Lambda Array')
        for j = 1
            if or(or(matches(data(i).type{j}, 'Full Extension'), matches(data(i).type{j}, 'Partial Extension')),matches(data(i).type{j}, 'Extension to 5pN')) 
                time_bool = and(data(i).time>data(i).times{j}(1), data(i).time<data(i).times{j}(2));
                plot(data(i).distance(time_bool), data(i).force(time_bool), 'Color', 'black', 'LineWidth', line_width)
                n_is_unwrap = n_is_unwrap+1;
            end
        end
        for j = 1:length(data(i).type)
            if or(matches(data(i).type{j}, 'Full Retraction'), matches(data(i).type{j}, 'Partial Retraction'))
                time_bool = and(data(i).time>data(i).times{j}(1), data(i).time<data(i).times{j}(2));
                plot(data(i).distance(time_bool), data(i).force(time_bool), 'Color', 'red', 'LineWidth', line_width)
                n_is_retrac = n_is_retrac+1;
            end
        end
    elseif  matches(data(i).condition, 'No Nucleosome Control')
       for j = 1:length(data(i).type)
            if matches(data(i).type{j}, 'Full Extension')
                time_bool = and(data(i).time>data(i).times{j}(1), data(i).time<data(i).times{j}(2));
                plot(data(i).distance(time_bool), data(i).force(time_bool), 'Color', 'green', 'LineWidth', line_width)
                n_is_refer = n_is_refer+1;
            elseif matches(data(i).type{j}, 'Full Retraction')
                time_bool = and(data(i).time>data(i).times{j}(1), data(i).time<data(i).times{j}(2));
                plot(data(i).distance(time_bool), data(i).force(time_bool), 'Color', 'green', 'LineWidth', line_width)
            end
        end 
    end
end

set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150])
xlim([9.9415   17.6577])
ylim([0   87.2909])

disp(newline)
disp(strcat('n = (unwrapping events): ', num2str(n_is_unwrap)))
disp(strcat('n = (retraction events): ', num2str(n_is_retrac)))
disp(strcat('n = (reference curves)', num2str(n_is_refer)))


%% Supplemental Figure 
% close all
figure('Name', 'Supplemental Figure')

% for i = 14:length(data)% 
i = 16;% take this one as an example 
bool = and(data(i).time<data(i).times{1}(2), data(i).time>data(i).times{1}(1));
force = data(i).force(bool);
distance = data(i).distance(bool);
time = data(i).time(bool);
plot(distance, force)
ylim([4.4435   20.5900])
xlim([13.2974   15.1299])
rectangle('Position',[14.6335 13 14.8220-14.6335 4])
hold on
xline(data(i).fivepNdis(1))
disp(i)
% pause
% end


%% (Supplemental Figure) 
figure('Name', 'Supplemental Figure')
% boolis = and(and(data(15).force<20.5, data(15).force>19.5), data(15).time <135);
boolis =and(force<15.5, force>14.5);
plot(time(boolis), distance(boolis))
% for i = 1:length(avgy)
%     yline(avgy(i))
% end
% zoom on
% pause
% xy = ginput(2);
xy = [73.4095   14.6649; 79.1513   14.6748];
boolis2 = and(time <xy(2,1), time>xy(1,1));
avgy = mean(distance(boolis2));
yline(avgy);
hold on
for i = 1:25
    yline(avgy+(0.025*(i-1)), 'Color', 'red')
end
ylim([ 14.6335   14.8220])
xlim([ 69.5161  108.8710])
copypaste_time = []; copypaste_distance=[];
copypaste_time = time(and(time>69.5161, time<108.8710))';
copypaste_distance = distance(and(distance>14.6335 , distance<14.8220))';
%% (Panel D) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except data summary
warning('off') 
figure('Name', 'Panel D')
hold on
counta = 1;
for i = [3, 7, 12]% this is for the last figure 
clearvars -except data summary i counta slopeval_holding yintercept_holding pvalue_holding rsquared_holding mean_MSD_from_hist x_save y_save
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
title('Figure 3: Panel X')
y_upperlim = 0.126258+0.0189;
ylim([0 y_upperlim])
xlim([0 2])
hold off
%% Supplemental Figure (Find Confined Area from mean MSD plot)
% close all
figure('Name', 'Supplementary Figure 2')

counta = 1;
for i = [3, 7, 12]%[2, 3, 1]
clearvars -except data summary i counta slopeval_holding yintercept_holding pvalue_holding rsquared_holding mean_MSD_from_hist x_save y_save
[M, F] = mode(summary(i).line_time);
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
% figure
% histogram(summary(i).dcoef(and(mode_restriction, fit_restriction)),30)
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

shadedErrorBar(xm, ym, ys)
xlabel('time lag (sec)')
ylabel('MSD (\mum^2)')
xlim([0 2])
hold on
% plot(mdl)
counta = counta+1;
% figure(7)
fit = fit(xm(xm<2)', ym(xm<2)', fittype('m*(1-exp(-b*x))+c'));
plot(fit, xm(xm<2), ym(xm<2), 'g') 
p12 = predint(fit,xm,0.95,'observation','on');
plot(xm,p12,'m--')
hold on
disp(strcat('asymptote is = ', num2str(abs(fit.m+fit.c)), ' meaning the confined length is: ', num2str(sqrt(abs(fit.m+fit.c)))))
disp(strcat(summary(i).condition, ' has an n = ', num2str(n)))
legend('off')
end
xlim([0 2])


%% Main figure Panel X

space = 3;
longest = 400;   %400 for SWR1 on nuc array
shortest = 100;   % 100 for SWR1 on nuc array
% close all
i = 7;
crop_all = summary(i).crop;
kymographs_all = summary(i).kymograph;
crop_area_all = summary(i).crop_coords;
line_time = summary(i).line_time;
slope = summary(i).slope;
pval = summary(i).pvalue;
r2 = summary(i).rsquared;
bool_keep = summary(i).bool_kept;
% Number of trajectories to plot:
num_2_plot = 15; 
% green map
green_map = zeros(3);
for i = 1:11
green_map(i,2) = i*0.1-0.1;
end
bw_map = zeros(3);
for i = 1:11
bw_map(i,1) = i*0.1-0.1;
bw_map(i,2) = i*0.1-0.1;
bw_map(i,3) = i*0.1-0.1;
end
length_crops = [];
width_crops = [];
crop = crop_all(bool_keep);
crop_area = crop_area_all(bool_keep);

kymo = kymographs_all(bool_keep);
lt = line_time(bool_keep);
for i = 1:length(crop)
length_crops(i,1) = size(crop{i},1);
width_crops(i,1)= size(crop{i}, 2);
end
[sorted,idx] = sort(length_crops,'descend');
s = find(sorted<longest & sorted >shortest);
w = 0;
for i = 1:num_2_plot
    j = idx(s(i));
    w = w + width_crops(j);
%     l = length_crops(j(1));
    line_time_plotted(i) = line_time(j);
end
mean_time = mean(line_time_plotted)/1000;% In seconds
j = idx(s(1));
l = length_crops(j);
z =[];
% PLOT WITH WHITE SPACE BETWEEN TRACES
z = ones(w+(num_2_plot*space)-space, l);
z_mask = zeros(w+(num_2_plot*space)-space, l);
position_for_plot = 0;
for i = 1:num_2_plot
    j = idx(s(i));
    if mean(mean(crop{j})) > 1
       for_plot = crop{j}./4;
    else 
       for_plot = crop{j};
    end
    z(position_for_plot+1:position_for_plot+width_crops(j),1:length_crops(j)) = rot90(for_plot);
    z_mask(position_for_plot+1:position_for_plot+width_crops(j),1:length_crops(j)) = ones(size(rot90(for_plot)));
    position_for_plot = position_for_plot + width_crops(j)+space;
    montage_cell{i} = rot90(for_plot);   
end
% figure
% imagesc(z_mask)
% x_ticks_are = xticks;
% new_xticks_are = round(x_ticks_are*mean_time);
% xticklabels(new_xticks_are)
% xlabel('time (s)')
% y_ticks_are = yticks;
% new_yticks_are = y_ticks_are*0.1;
% ylabel('distance (\mum)')
% yticklabels(new_yticks_are)
% ax = gca;
% ax.FontSize = 16;
% axis square
% colormap(bw_map)
% caxis([0 4])
figure
imagesc(z)
x_ticks_are = xticks;
new_xticks_are = round(x_ticks_are*mean_time);
xticklabels(new_xticks_are)
xlabel('time (s)')
y_ticks_are = yticks;
new_yticks_are = y_ticks_are*0.1;
ylabel('distance (\mum)')
yticklabels(new_yticks_are)
ax = gca;
ax.FontSize = 16;
axis square
colormap(green_map)
caxis([0 4])

