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

%% Convert to a table
% nrow=length(colonies);
% times=nan(nrow, min);
% lengths=nan(nrow, min);
% intensities=nan(nrow, min);
% norm_intensities=nan(nrow, min);
% A=nan(nrow, 1);
% B=nan(nrow, 1);
% C=nan(nrow, 1);
% D=nan(nrow, 1);

times=[];
lengths=[];
intensities=[];
norm_intensities=[];
A=[];
B=[];
C=[];
D=[];

for i=1:length(data)
    i
    
    time=repmat(data(i).time(1:min), height(data(i).length),1);
    times=[times; time];
    
    lengths=[lengths; data(i).length(:,1:min)];
    intensities=[intensities; data(i).intensity(:, 1:min)];
    norm_intensities=[norm_intensities; data(i).norm_intensity(:, 1:min)];
    A=[A; data(i).a(:, 1)];
    B=[B; data(i).b(:, 1)];
    C=[C; data(i).c(:, 1)];
    D=[D; data(i).d(:, 1)];   
end

dataTable=table(experiments, colonies, times, lengths, intensities, norm_intensities, A, B, C, D, 'VariableNames', {'experiment', 'colony', 'time', 'length', 'intensity', 'norm_intensity','a', 'b', 'c', 'd'});
%% plot data
% figure(1), hold on
% for i=1:height(data)
%     if strcmp(data(i).experiment, '07242021_Exp1')==1
%         plot(data(i).time(1:min), data(i).norm_intensity(:, 1:min), '-m')
%     elseif strcmp(data(i).experiment, '08172021_Exp1')==1
%         plot(data(i).time(1:min), data(i).norm_intensity(:, 1:min), '-c')
%     end
% end
% 
% figure(2), hold on
% for i=1:height(data)
%     if strcmp(data(i).experiment, '08172021_Exp1')==1
%         plot(data(i).time(1:min), data(i).norm_intensity(:, 1:min), '-c')
%     end
% end

%% save data
cd(dirsave)
save('fitting.mat')
%writetable(dataTable, 'dataTable.csv')

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