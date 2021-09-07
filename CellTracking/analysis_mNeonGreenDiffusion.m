%%PlotScrapBook.m
%Purpose: to save code for one-off plots I've made

clear, close all

%% user input
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Downloads/NYU/Year3_2021_Fall/presentations/09062021";

data = struct('experiment', [], 'colony', struct('norm_green', [], 'time', [], 'a', [], 'b', [], 'c', [], 'd', [], 'fitModel', []));
cd(datadir)
cutoff=13;

%% LB, frame rate = 2
data(1).experiment = '07162021_Exp3';
load("07162021_Exp3_colony2_dm.mat", 'dataTable', 'time')
norm_green = table2array(dataTable(dataTable.halfie==0, 3));
%data(1).colony(1).norm_green = norm_green(:, cutoff:end);
data(1).colony(1).norm_green = norm_green;
time = time-time(1);
%data(1).colony(1).time = time(cutoff:end);
data(1).colony(1).time = time;

load("07162021_Exp3_colony3_dm.mat", 'dataTable')
norm_green = table2array(dataTable(dataTable.halfie==0, 3));
%data(1).colony(2).norm_green = norm_green(:, cutoff:end);
data(1).colony(2).norm_green = norm_green;
%data(1).colony(2).time = time(cutoff:end);
data(1).colony(2).time = time;

%% LB, frame rate = 1
data(2).experiment = '08262021_Exp1';
load("08262021_Exp1_colony1_dm.mat", 'dataTable', 'time')
data(2).colony(1).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(2).colony(1).time = time;

load("08262021_Exp1_colony2_dm.mat", 'dataTable', 'time')
data(2).colony(2).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(2).colony(2).time = time;
 
load("08262021_Exp1_colony3_dm.mat", 'dataTable', 'time')
data(2).colony(3).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(2).colony(3).time = time;

load("08262021_Exp1_colony4_dm.mat", 'dataTable', 'time')
data(2).colony(4).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(2).colony(4).time = time;

load("08262021_Exp1_colony5_dm.mat", 'dataTable', 'time')
data(2).colony(5).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(2).colony(5).time = time;

load("08262021_Exp1_colony6_dm.mat", 'dataTable', 'time')
data(2).colony(6).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(2).colony(6).time = time;

%% EDTA
data(3).experiment = '07162021_Exp4';
load("07162021_Exp4_colony1_dm.mat", 'dataTable', 'time')
data(3).colony(1).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(3).colony(1).time = time-time(1);

% load("07162021_Exp4_colony2_dm.mat", 'dataTable', 'time')
% data(3).colony(2).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
% data(3).colony(2).time = time;

load("07162021_Exp4_colony3_dm.mat", 'dataTable', 'time')
data(3).colony(3).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(3).colony(3).time = time-time(1);

load("07162021_Exp4_colony4_dm.mat", 'dataTable', 'time')
data(3).colony(4).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(3).colony(4).time = time-time(1);

%% Mg2+
data(4).experiment = '07152021_Exp1';
load("07152021_Exp1_colony1_dm.mat", 'dataTable', 'time')
data(4).colony(1).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(4).colony(1).time = time-time(1);

load("07152021_Exp1_colony2_dm.mat", 'dataTable', 'time')
data(4).colony(2).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(4).colony(2).time = time-time(1);

load("07152021_Exp1_colony3_dm.mat", 'dataTable', 'time')
data(4).colony(3).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(4).colony(3).time = time-time(1);
%% PBS, frame rate = 1
data(5).experiment = '08272021_Exp1';

load("08272021_Exp1_colony1_dm.mat", 'dataTable', 'time')
data(5).colony(1).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(5).colony(1).time = time;

load("08272021_Exp1_colony2_dm.mat", 'dataTable', 'time')
data(5).colony(2).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(5).colony(2).time = time;

load("08272021_Exp1_colony3_dm.mat", 'dataTable', 'time')
data(5).colony(3).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(5).colony(3).time = time;

load("08272021_Exp1_colony4_dm.mat", 'dataTable', 'time')
data(5).colony(4).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(5).colony(4).time = time;

load("08272021_Exp1_colony5_dm.mat", 'dataTable', 'time')
data(5).colony(5).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(5).colony(5).time = time;

load("08272021_Exp1_colony6_dm.mat", 'dataTable', 'time')
data(5).colony(6).norm_green = table2array(dataTable(dataTable.halfie==0, 3));
data(5).colony(6).time = time;

%% fit curves
for i=1:length(data)
    for j=1:length(data(i).colony)
        [data(i).colony(j).a, data(i).colony(j).b, data(i).colony(j).c, data(i).colony(j).d, data(i).colony(j).fitModel] = expFitting(data(i).colony(j).norm_green, data(i).colony(j).time);
    end
end

%% plot data
xdata = [0:180];
for i=1:length(data)
    for j=1:length(data(i).colony)
%             figure(1), hold on
%             histogram(data(i).colony(j).a),
%             histogram(data(i).colony(j).c)
%             legend({'a', 'c'})
%             hold off 
%             
%             figure(2), hold on
%             histogram(data(i).colony(j).b), 
%             histogram(data(i).colony(j).d)
%             legend({'b', 'd'})
%             hold off
            
%             figure(3), hold on
%             for k=1:length(data(i).colony(j).fitModel)
%                 if isempty(data(i).colony(j).norm_green(k,:))==0 & isnan(data(i).colony(j).norm_green(k,:))==0
%                     plot(data(i).colony(j).fitModel{k,1}, data(i).colony(j).time, data(i).colony(j).norm_green(k,:));
%                     title([data(i).experiment, ', colony ', num2str(j)])
%                     legend('off')
%                 end
%             end
%             pause, close all
%             
              figure(4), hold on
              for k=1:length(data(i).colony(j).fitModel)
                if isnan(data(i).colony(j).norm_green(k,:))==0
                  ydata = data(i).colony(j).fitModel{k,1}(xdata);
                if i==1
                  %plot(xdata, ydata, '-g')
                elseif i==2
                  plot(xdata, ydata, '-r')
                elseif i==3
                    %plot(xdata, ydata, '-m')
                elseif i==4
                    %plot(xdata, ydata, '-c')
                elseif i==5
                    plot(xdata, ydata, '-b')
                end
                
              %title('Intensity vs Time, red=LB, blue=LB2, magenta=EDTA, cyan=Mg2+, green=PBS')
              title('Intensity vs Time, red=LB, blue=PBS')
              %pause, close all    
                end
              end
              
    end
end
%% save data
cd(dirsave)
save('fitting.mat')

%% individual experiment plots
for i=1 %1:length(data)
    for j=1:length(data(i).colony)
            figure(1), hold on
            for k=1:length(data(i).colony(j).fitModel)
                if isempty(data(i).colony(j).norm_green(k,:))==0 & isnan(data(i).colony(j).norm_green(k,:))==0
                    %plot(data(i).colony(j).fitModel{k,1}, data(i).colony(j).time, data(i).colony(j).norm_green(k,:));
                    plot(data(i).colony(j).time, data(i).colony(j).norm_green(k,:));
                    title([data(i).experiment, ', colony ', num2str(j)])
                    legend('off')
                end
            end
            %pause, close all
    end
end
%% functions
% [a, b, c, d, fitScore]=expFitting(data(i).colony(j).norm_green, data(i).colony(j).time);

function [a, b, c, d, fitModel] = expFitting(norm_green, time)

    fitModel=cell(height(norm_green), 1);
    a=nan(height(norm_green), 1);
    b=nan(height(norm_green), 1);
    c=nan(height(norm_green), 1);
    d=nan(height(norm_green), 1);
    %fitScore=nan(height(norm_green), 1);
    
    for n=1:height(norm_green)
        if isnan(norm_green(n, 1))==0
            [xData, yData]=prepareCurveData(time, norm_green(n,:));
            %figure
            fitModel{n,1}=fit(xData, yData, 'exp2'); %fit to a linear exponential
            
            %if isnan(fitModel{n,1})==0    
                a(n,1) = fitModel{n,1}.a;
                b(n,1) = fitModel{n,1}.b;
                c(n,1) = fitModel{n,1}.c;
                d(n,1) = fitModel{n,1}.d;
    %             plot(f{n,1}, xData, yData)
    %             prompt = 'Q1: Is this a good fit? 0=No, 1=Yes ';
    %             fitScore(n,1) = input(prompt);
    %             close
            %else
             %   fitModel{n,1}=NaN;
            %end
        else
            fitModel{n,1}=NaN;
            %fitScore(n,1)=NaN;
            a(n,1)=NaN;
            b(n,1)=NaN;
            c(n,1)=NaN;
            d(n,1)=NaN;    
        end
    end
end