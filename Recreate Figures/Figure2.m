%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the following datasets:
% All Diffusion Dataset (variable name: summary)

%% (Panel C) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
close all
figure('Name', 'Panel C')
clearvars -except summary
kymo = imread('C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\03 My Publications\01_SWR1_sliding\Data Per Figure\MF2 Panel C\20190829-201214 Kymograph 10(contrast corrected).tif');
% kymo = imread('C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\03 My Publications\01_SWR1_sliding\Data Per Figure\MF2 Panel C\20190829-201214 Kymograph 10 (before contrast correction).tiff');
kymo_ID = size(kymo, 1)*size(kymo, 2);  
for i = 14
    for j = 1:length(summary(i).kymograph)
        summary(i).ID(j) = size(summary(i).kymograph{j}, 1)*size(summary(i).kymograph{j}, 2);        
    end
end
for i = 14
    for j = 1:length(summary(i).ID)
        if summary(i).ID(j) == kymo_ID
            disp(j);
            disp(summary(i).line_time(j));
            kymo_LT = summary(i).line_time(j);
        end
    end
end          
imagesc(kymo)
old = [];
old = get(gca,'XTick');
new = [];
new = round(old.*kymo_LT/1000);
xticklabels(new)
old_y = [];
old_y = get(gca,'YTick');
new_y = [];
new_y = flip(old_y)/10;
yticklabels(new_y)
xlabel('time(s)')
ylabel('position(um)')
title(['RGB kymograph: ID # ', num2str(kymo_ID)])
xlim([0 487])
    
%% Figure 2 (Panel D) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
i = 14; % no ATP example  
len = cellfun(@length, summary(i).MSDs);
figure('Name', 'Panel D')
subplot(1, 2, 1)
range_search = and(len>100, len<200);
range_search = and(range_search, summary(i).bool_kept); % apply fitting criteria
range_search = and(range_search, summary(i).outlier_exclusion); % apply outlier criteria 

% get a representative individual MSD from within each diffusion
% coefficient regime, and of a particular length

low1 = summary(14).dcoef<0.01;
low2 = and(summary(14).dcoef>0.01, summary(14).dcoef<0.02);
low3 = and(summary(14).dcoef>0.02, summary(14).dcoef<0.03);
low4 = and(summary(14).dcoef>0.03, summary(14).dcoef<0.04);
low5 = and(summary(14).dcoef>0.04, summary(14).dcoef<0.05);

ex1 = and(range_search, low1);
ex2 = and(range_search, low2);
ex3 = and(range_search, low3);
ex4 = and(range_search, low4);
ex5 = and(range_search, low5);

aex1 = and(ex1, summary(14).bool_kept);
aex2 = and(ex2, summary(14).bool_kept);
aex3 = and(ex3, summary(14).bool_kept);
aex4 = and(ex4, summary(14).bool_kept);
aex5 = and(ex5, summary(14).bool_kept);
aex5(find(aex5, 1, 'first')) =0 ;

x = [0:10];
hold on

plot(summary(14).MSDs{find(aex1, 1, 'first')}(:,1),summary(14).MSDs{find(aex1, 1, 'first')}(:,2),'Color', '#0072BD')
y = summary(14).slope(find(aex1, 1, 'first'))*x+summary(14).yintercept(find(aex1, 1, 'first'));
plot(x,y,'Color', '#0072BD')
copypastethis_time = [];
copypastethis_MSD = [];
disp(strcat('m = ', num2str(summary(14).slope(find(aex1, 1, 'first')))))
disp(strcat('b = ', num2str(summary(14).yintercept(find(aex1, 1, 'first')))))
disp(strcat('r2 = ', num2str(summary(14).rsquared(find(aex1, 1, 'first')))))
disp(strcat('pval = ', num2str(summary(14).pvalue(find(aex1, 1, 'first')))))
copypastethis_time = summary(14).MSDs{find(aex1, 1, 'first')}(:,1);
copypastethis_MSD = summary(14).MSDs{find(aex1, 1, 'first')}(:,2);

plot(summary(14).MSDs{find(aex2, 1, 'first')}(:,1),summary(14).MSDs{find(aex2, 1, 'first')}(:,2),'Color', '#D95319')
y = summary(14).slope(find(aex2, 1, 'first'))*x+summary(14).yintercept(find(aex2, 1, 'first'));
plot(x,y,'Color', '#D95319')
disp(newline)
copypastethis_time = [];
copypastethis_MSD = [];
disp(strcat('m = ', num2str(summary(14).slope(find(aex2, 1, 'first')))))
disp(strcat('b = ', num2str(summary(14).yintercept(find(aex2, 1, 'first')))))
disp(strcat('r2 = ', num2str(summary(14).rsquared(find(aex2, 1, 'first')))))
disp(strcat('pval = ', num2str(summary(14).pvalue(find(aex2, 1, 'first')))))
copypastethis_time = summary(14).MSDs{find(aex2, 1, 'first')}(:,1);
copypastethis_MSD = summary(14).MSDs{find(aex2, 1, 'first')}(:,2);

plot(summary(14).MSDs{find(aex3, 1, 'first')}(:,1),summary(14).MSDs{find(aex3, 1, 'first')}(:,2),'Color', '#EDB120')
y = summary(14).slope(find(aex3, 1, 'first'))*x+summary(14).yintercept(find(aex3, 1, 'first'));
plot(x,y,'Color', '#EDB120')
disp(newline)
copypastethis_time = [];
copypastethis_MSD = [];
disp(strcat('m = ', num2str(summary(14).slope(find(aex3, 1, 'first')))))
disp(strcat('b = ', num2str(summary(14).yintercept(find(aex3, 1, 'first')))))
disp(strcat('r2 = ', num2str(summary(14).rsquared(find(aex3, 1, 'first')))))
disp(strcat('pval = ', num2str(summary(14).pvalue(find(aex3, 1, 'first')))))
copypastethis_time = summary(14).MSDs{find(aex3, 1, 'first')}(:,1);
copypastethis_MSD = summary(14).MSDs{find(aex3, 1, 'first')}(:,2);

plot(summary(14).MSDs{find(aex4, 1, 'first')}(:,1),summary(14).MSDs{find(aex4, 1, 'first')}(:,2),'Color', '#7E2F8E')
y = summary(14).slope(find(aex4, 1, 'first'))*x+summary(14).yintercept(find(aex4, 1, 'first'));
plot(x,y,'Color', '#7E2F8E')
disp(newline)
copypastethis_time = [];
copypastethis_MSD = [];
disp(strcat('m = ', num2str(summary(14).slope(find(aex4, 1, 'first')))))
disp(strcat('b = ', num2str(summary(14).yintercept(find(aex4, 1, 'first')))))
disp(strcat('r2 = ', num2str(summary(14).rsquared(find(aex4, 1, 'first')))))
disp(strcat('pval = ', num2str(summary(14).pvalue(find(aex4, 1, 'first')))))
copypastethis_time = summary(14).MSDs{find(aex4, 1, 'first')}(:,1);
copypastethis_MSD = summary(14).MSDs{find(aex4, 1, 'first')}(:,2);

plot(summary(14).MSDs{find(aex5, 1, 'first')}(:,1),summary(14).MSDs{find(aex5, 1, 'first')}(:,2),'Color', '#77AC30')
y = summary(14).slope(find(aex5, 1, 'first'))*x+summary(14).yintercept(find(aex5, 1, 'first'));
plot(x,y,'Color', '#77AC30')
disp(newline)
copypastethis_time = [];
copypastethis_MSD = [];
disp(strcat('m = ', num2str(summary(14).slope(find(aex5, 1, 'first')))))
disp(strcat('b = ', num2str(summary(14).yintercept(find(aex5, 1, 'first')))))
disp(strcat('r2 = ', num2str(summary(14).rsquared(find(aex5, 1, 'first')))))
disp(strcat('pval = ', num2str(summary(14).pvalue(find(aex5, 1, 'first')))))
copypastethis_time = summary(14).MSDs{find(aex5, 1, 'first')}(:,1);
copypastethis_MSD = summary(14).MSDs{find(aex5, 1, 'first')}(:,2);

hold off
xlim([0 8])
ylim([0 .8])
xlabel('time(s)')
ylabel('MSD(um2/sec)')

% Plot only 2 seconds 
subplot(1, 2, 2) 
x = [0:10];
hold on
plot(summary(14).MSDs{find(aex1, 1, 'first')}(:,1),summary(14).MSDs{find(aex1, 1, 'first')}(:,2),'Color', '#0072BD')
y = summary(14).slope(find(aex1, 1, 'first'))*x+summary(14).yintercept(find(aex1, 1, 'first'));
plot(x,y,'Color', '#0072BD')

plot(summary(14).MSDs{find(aex2, 1, 'first')}(:,1),summary(14).MSDs{find(aex2, 1, 'first')}(:,2),'Color', '#D95319')
y = summary(14).slope(find(aex2, 1, 'first'))*x+summary(14).yintercept(find(aex2, 1, 'first'));
plot(x,y,'Color', '#D95319')

plot(summary(14).MSDs{find(aex3, 1, 'first')}(:,1),summary(14).MSDs{find(aex3, 1, 'first')}(:,2),'Color', '#EDB120')
y = summary(14).slope(find(aex3, 1, 'first'))*x+summary(14).yintercept(find(aex3, 1, 'first'));
plot(x,y,'Color', '#EDB120')

plot(summary(14).MSDs{find(aex4, 1, 'first')}(:,1),summary(14).MSDs{find(aex4, 1, 'first')}(:,2),'Color', '#7E2F8E')
y = summary(14).slope(find(aex4, 1, 'first'))*x+summary(14).yintercept(find(aex4, 1, 'first'));
plot(x,y,'Color', '#7E2F8E')

plot(summary(14).MSDs{find(aex5, 1, 'first')}(:,1),summary(14).MSDs{find(aex5, 1, 'first')}(:,2),'Color', '#77AC30')
y = summary(14).slope(find(aex5, 1, 'first'))*x+summary(14).yintercept(find(aex5, 1, 'first'));
plot(x,y,'Color', '#77AC30')

hold off
xlim([0 2])
ylim([0 .2])
xlabel('time(s)')
ylabel('MSD(um2/sec)')

sgtitle('Panel D')
set(gcf, 'Position', [373.0000  284.3333  796.3333  341.0000])
%% Figure 2 (Panel E) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Panel E')

for i = [12, 14]
    if i == 12
        subplot(2, 1, 1)
    elseif i == 14
        subplot(2, 1, 2)
    end
criteria = and(summary(i).bool_kept,  summary(i).outlier_exclusion);
copypaste = [];
histogram(summary(i).dcoef(criteria), 20, 'BinLimits', [0 0.14], 'Normalization', 'probability', 'DisplayStyle', 'bar', 'DisplayName', strcat(summary(i).condition, ' n = ', num2str(length(summary(i).dcoef(criteria)))))
copypaste = summary(i).dcoef(criteria)';
legend
title(summary(i).condition)
xlabel("Dcoff \mum^2/sec")
ylabel('Probability')
end
set(gcf, 'Position', [360    52   308   566])
%% (Panel F) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% no ATP representative %%%%%%%%%%%%%%%%%%%%%%%%
clc
space = 3;
longest = 250;   %400 for SWR1 on nuc array
shortest = 100;   % 100 for SWR1 on nuc array
% close all
i = 14;
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




%%%%%%%%%%%%%%%%%%%%% dCas9 representative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
space = 3
longest = 170   %400 for SWR1 on nuc array
shortest = 10   % 100 for SWR1 on nuc array
% close all
i = 12;
crop_all = summary(i).crop;
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
% black and white map
bw_map = zeros(3);
for i = 1:11
bw_map(i,1) = i*0.1-0.1;
bw_map(i,2) = i*0.1-0.1;
bw_map(i,3) = i*0.1-0.1;
end
% pvalL = pval < 0.1;
% r2L = r2 > 0.8;
length_crops = [];
width_crops = [];
crop = crop_all(bool_keep);
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
z = ones(w+(num_2_plot*space)-space, l);
z_mask = zeros(w+(num_2_plot*space)-space, l);
position_for_plot = 0;
for i = 1:num_2_plot
    j = idx(s(i));   
    if mean(mean(crop{j})) > 5
       for_plot = crop{j}./8;
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
caxis([0 8])
