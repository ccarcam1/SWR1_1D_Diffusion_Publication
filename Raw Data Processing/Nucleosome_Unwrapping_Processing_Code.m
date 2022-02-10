%% DATA PROCESSING, BEFORE PLOTS
% Step 1: Sort Force Distance Plots into extensions and retraction
close all
clc
hold off
for i = 1:length(data)
    disp(i)
    plot(data(i).time, data(i).force, "LineWidth", 0.25, 'Color', 'black')
    xlabel('Time')
    ylabel('Force')
    count = 1;
    keep_going = 1;
    while keep_going == 1
    pause
    answer = questdlg('Extension or Retraction', ...
	'Question 1', ...
	'Full Extension','Full Retraction','Neither', 'Neither')
    % Handle response
        switch answer
            case {'Full Extension','Full Retraction'}
                disp([answer])
                x = []; 
                y = [];
                [x,y] = ginput(2);
                data(i).times{count} = x;
                data(i).type{count} = answer;
            case 'Neither'
                answer2 = questdlg('Extension or Retraction', ...
                'Question 2', ...
                'Partial Extension','Partial Retraction', 'Neither','Neither')
                switch answer2 
                    case {'Partial Extension','Partial Retraction'}
                    disp([answer2])
                    x = []; 
                    y = [];
                    [x,y] = ginput(2);
                    data(i).times{count} = x;
                    data(i).type{count} = answer2;
                     case 'Neither'    
                     answer3 = questdlg('Extension or Retraction', ...
                    'Question 2', ...
                    'Extension to 5pN','Retraction from 5pN','Bad Curve','Bad Curve')
                        switch answer3
                            case {'Extension to 5pN','Retraction to 5pN'}
                                disp([answer3])
                                x = []; 
                                y = [];
                                [x,y] = ginput(2);
                                data(i).times{count} = x;
                                data(i).type{count} = answer3;
                            case 'Bad Curve' 
                                data(i).times{count} = [];
                                data(i).type{count} = answer3; 
                                break
                        end
                end
        end
        answer4 = questdlg('Keep Going?', ...
        'Question 2', ...
        'Keep Going','Next Up','Next Up')
 % Handle response
        switch answer4
            case 'Keep Going'
                keep_going = 1;
            case 'Next Up'
                keep_going = 0;
        end
    count = count +1;    
    end
end

%% Step 2 Count the number of nucleosome unwrapping events by eye 
close all
which_used = [];
k = 1;
for i = 15:length(data)
    if not(matches(data(i).type{1}, 'Extension to 5pN'))
    plot(data(i).distance(and(data(i).time > data(i).times{1}(1), data(i).time < data(i).times{1}(2))), data(i).force(and(data(i).time > data(i).times{1}(1), data(i).time < data(i).times{1}(2))), "LineWidth", 0.25, 'Color', 'black')
    pause
    answer = questdlg('Good or Bad', ...
        'Question 1', ...
        'Good','Bad','Neither', 'Neither')
        % Handle response
        switch answer
            case "Good"
                unwrapping(k) = input("How many unwrapping events: ");
                data(i).useunwrappingcounts2 = 1;
                data(i).unwrappingcounts2 = unwrapping(k);
            case "Bad"
                data(i).useunwrappingcounts2 = 0;
        end

        which_used(k) = i;
        k = k+1;
        disp(i)
%         pause
    end 
end
hold off
xlim([10 17.5])
ylim([0 65])


%% Get the 5pN distance and find compaction ratio 
% first determine computationally
for i = 25:length(data)
    for j = 1:length(data(i).type)
        plot(data(i).distance(and(data(i).time > data(i).times{j}(1), data(i).time < data(i).times{j}(2))), data(i).force(and(data(i).time > data(i).times{j}(1), data(i).time < data(i).times{j}(2))), "LineWidth", 0.25, 'Color', 'black')
        force = [];
        force = data(i).force(and(data(i).time > data(i).times{j}(1), data(i).time < data(i).times{j}(2)));
        distance = [];
        distance = data(i).distance(and(data(i).time > data(i).times{j}(1), data(i).time < data(i).times{j}(2)));
        data(i).fivepNdis(j) = mean(distance(round(force) == 5));
        xline(data(i).fivepNdis(j))
        zoom on
        prompt = 'Does the 5pN look good?';
        answeris = input(prompt, 's')
        if not(matches(answeris, 'Yes'))
           xy =  ginput(1);
           data(i).fivepNdis(j) = xy(1);
           xline(data(i).fivepNdis(j), 'Color', 'red')
           pause
        end
    end
end