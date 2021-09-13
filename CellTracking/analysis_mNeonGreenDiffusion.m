%%PlotScrapBook.m
%Purpose: to save code for one-off plots I've made

clear, close all

%% user input
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09142021_analysis";

data = struct('experiment', [], 'colony', [], 'length', [], 'intensity', [], 'norm_intensity', [], 'time', [], 'a', [], 'b', [], 'c', [], 'd', [], 'fitModel', []);
cd(datadir)

files=dir('*.mat');

experiments=[];
colonies=[];

min=6000;
for i=1:height(files)
    
    data(i).experiment=files(i).name(1:end-15);
    data(i).colony=files(i).name(15:21);
    load(files(i).name, 'dataTable', 'time')
    data(i).length=table2array(dataTable(dataTable.halfie==0, 1));
    data(i).intensity=table2array(dataTable(dataTable.halfie==0, 2));
    data(i).norm_intensity=table2array(dataTable(dataTable.halfie==0, 3));
    
    if time(1)~=0 
        time = time-time(1);
    end
  
    data(i).time=time;
    
    if length(time) < min
        min=length(time);
    end
    
    line1=convertCharsToStrings(data(i).experiment);
    line2=convertCharsToStrings(data(i).colony);
    
    mat1=repmat(line1, height(data(i).length),1);
    mat2=repmat(line2, height(data(i).length),1);
    
    experiments=[experiments; mat1];
    colonies=[colonies;mat2];
    
end

%% Convert to a table

%% Notes
%Mg2+, time = [0:2:180]
%LB, frame rate = 2, time = [0:2:180]
%EDTA, time = [0:2:180]
%LB + sodium azide, frame rate = 2, time = [0:2:180]
%LB + sodium azide, frame rate = 2, time = [0:2:180]
%LB, frame rate = 1, time = [0:91]
%PBS, frame rate = 1, incubation = 2 min, time = [0:90]
%PBS, frame rate = 1, incubation = 8 min, time = [0:1:45]
%PBS, frame rate = 1, incubation = 16 min, time = [0:1:45]

%% fit data to exp2
for i=1:length(data)
    [data(i).a, data(i).b, data(i).c, data(i).d, data(i).fitModel] = expFitting(data(i).norm_intensity, data(i).time);
end

%% plot data
% figure(1), hold on
% for i=16:21
%     plot(data(i).time(1:min), data(i).norm_green(:,1:min), '-r')
% end
% 
% for i=22:27
%     plot(data(i).time(1:min), data(i).norm_green(:,1:min), '-b')
% end
% 
% for i=28:32
%     plot(data(i).time(1:min), data(i).norm_green(:,1:min), '-c')
% end

% figure(2), hold on
% for i=10:15
% %     for j=1:height(data(i).norm_green)
% %       plot(data(i).fitModel{j}, data(i).time, data(i).norm_green(j,:))
% %     end
% 
%     histogram(data(i).a)
% end


%% save data
cd(dirsave)
save('fitting.mat')

%% functions
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