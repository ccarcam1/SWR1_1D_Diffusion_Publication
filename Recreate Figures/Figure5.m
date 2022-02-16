%% Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the following datasets:
% All Diffusion Dataset (variable name: summary)
% dCas9 Bypass Dataset (variable name: corr)

%% Reproduce the pie-chart

clearvars -except summary corr 
clc
close all
figure('Name', 'Panel C')

% The message section is where I separate colocalized traces into either:
% 1: Stuck (Immobile)
% 2: Stuck + Bounce (Some immobile segments with confined diffusion)
% 3: Bounce (Confined diffusion) 
% 4: Cross over: After bleach
% 5: Cross over: Before bleach

message = {corr.message}.';
counter = 1;
pie_data = [];
for i = 1:length(message)
    if ~isempty(message{i})
        for j = 1:size(message{i},1)
            pie_data(counter) = message{i}{j} ;
            counter = counter + 1;
        end
    end
end
 
% Create pie charts
[GC,GR] = groupcounts(pie_data');
GCC = [GC(1) GC(2)+GC(3) GC(4)+GC(5)];
% length_list = length(pie_data);
% percentages = (GC./length_list).*100;
pie(GCC)
% labels = {"Stuck", "Stuck + Bounce", "Bounce", "Cross over: After bleach", "Cross over: Before bleach"};
labels = {"Stuck", "Bounce", "Bypass"};
% Create legend
lgd = legend(labels);
disp(strcat('n = ', num2str(length(pie_data))))