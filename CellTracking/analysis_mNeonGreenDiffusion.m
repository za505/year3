%%PlotScrapBook.m
%Purpose: to save code for one-off plots I've made

clear, close all

%% user input
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09142021_analysis";
presdir="/Users/zarina/Downloads/NYU/Year3_2021_Fall/updates/09142021";

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

%annotation=["Mg2+", "LB, fr=2", "EDTA", "sodium azide", "sodium azide", "LB, fr=1", "PBS, 2 min", "PBS, 8 min", "PBS, 16 min"];
SA1=116:123;
SA2=124:151;
PBS1=239:340;
PBS2=341:417;
PBS3=418:462;
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

dataTable=table(categorical(experiments), categorical(colonies), times, lengths, intensities, norm_intensities, A, B, C, D, 'VariableNames', {'experiment', 'colony', 'time', 'length', 'intensity', 'norm_intensity','a', 'b', 'c', 'd'});
%% plot sodium azide data
cd(presdir)

figure(1), hold on
for i=SA1
    plot(times(i,:), lengths(i,:), '-r')
end
for i=SA2
    plot(times(i,:), lengths(i,:), '-b')
end
xlabel('Time (minutes)')
ylabel('Length (\mum)')
title('Length vs Time, Sodium Azide (07242021 in red, 08172021 in blue)')
saveas(gcf, 'length_SA.png')
saveas(gcf, 'length_SA.fig')

figure(2), hold on
for i=SA1
    plot(times(i,:), norm_intensities(i,:), '-r')
end
for i=SA2
    plot(times(i,:), norm_intensities(i,:), '-b')
end
xlabel('Time (minutes)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time, Sodium Azide (07242021 in red, 08172021 in blue)')
saveas(gcf, 'normintensity_SA.png')
saveas(gcf, 'normintensity_SA.fig')

figure(3), hold on
for i=SA1
    plot(times(i,:), intensities(i,:), '-r')
end
for i=SA2
    plot(times(i,:), intensities(i,:), '-b')
end
xlabel('Time (minutes)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time, Sodium Azide (07242021 in red, 08172021 in blue)')
saveas(gcf, 'intensity_SA.png')
saveas(gcf, 'intensity_SA.fig')

close all
%% plot PBS data
cd(presdir)

figure(1), hold on
for i=PBS1
    plot(times(i,:), lengths(i,:), '-r')
end
for i=PBS2
    plot(times(i,:), lengths(i,:), '-b')
end
for i=PBS3
    plot(times(i,:), lengths(i,:), '-c')
end
xlabel('Time (minutes)')
ylabel('Length (\mum)')
title('Length vs Time, 2 min=red, 8 min=blue, 16 min=cyan')
saveas(gcf, 'length_PBS.png')
saveas(gcf, 'length_PBS.fig')

figure(2), hold on
for i=PBS1
    plot(times(i,:), norm_intensities(i,:), '-r')
end
for i=PBS2
    plot(times(i,:), norm_intensities(i,:), '-b')
end
for i=PBS3
    plot(times(i,:), norm_intensities(i,:), '-c')
end
xlabel('Time (minutes)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time, 2 min=red, 8 min=blue, 16 min=cyan')
saveas(gcf, 'normintensity_PBS.png')
saveas(gcf, 'normintensity_PBS.fig')

figure(3), hold on
for i=PBS1
    plot(times(i,:), intensities(i,:), '-r')
end
for i=PBS2
    plot(times(i,:), intensities(i,:), '-b')
end
for i=PBS3
    plot(times(i,:), intensities(i,:), '-c')
end
xlabel('Time (minutes)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time, 2 min=red, 8 min=blue, 16 min=cyan')
saveas(gcf, 'intensity_PBS.png')
saveas(gcf, 'intensity_PBS.fig')

%% Histogram plots of PBS alpha
flagged=find(A>1.5|A<0);
figure(4), hold on
for i=PBS1
    if ismember(i, flagged)==0
        h1=histogram(A(i,:))
        h1.FaceColor = 'r';
        h1.EdgeColor = 'k';
        h1.NumBins = 1;
        h1.BinWidth = 0.1;
    end
end
xlabel('alpha')
ylabel('Counts')
title('Distribution of alpha coefficients, 2 min=red, 8 min=blue, 16 min=cyan')
saveas(gcf, 'hist_A_PBS1.png')
saveas(gcf, 'hist_A_PBS1.fig')

figure(5), hold on
for i=PBS2
    if ismember(i, flagged)==0
        h2=histogram(A(i,:))
        h2.FaceColor = 'b';
        h2.EdgeColor = 'k';
        h2.NumBins = 1;
        h2.BinWidth = 0.1;
    end
end
xlabel('alpha')
ylabel('Counts')
title('Distribution of alpha coefficients, 2 min=red, 8 min=blue, 16 min=cyan')
saveas(gcf, 'hist_A_PBS2.png')
saveas(gcf, 'hist_A_PBS2.fig')

figure(6), hold on
for i=PBS3
    if ismember(i, flagged)==0
        h3=histogram(A(i,:))
        h3.FaceColor = 'c';
        h3.EdgeColor = 'k';
        h3.NumBins = 1;
        h3.BinWidth = 0.1;
    end
end
xlabel('alpha')
ylabel('Counts')
title('Distribution of alpha coefficients, 2 min=red, 8 min=blue, 16 min=cyan')
saveas(gcf, 'hist_A_PBS3.png')
saveas(gcf, 'hist_A_PBS3.fig')

close all
%% Histogram plots of PBS beta
flagged=find(A>1.5|A<0);
figure(4), hold on
for i=PBS1
    if ismember(i, flagged)==0
        h1=histogram(B(i,:))
        h1.FaceColor = 'r';
        h1.EdgeColor = 'k';
        h1.NumBins = 1;
        h1.BinWidth = 0.1;
    end
end

for i=PBS2
    if ismember(i, flagged)==0
        h2=histogram(B(i,:))
        h2.FaceColor = 'b';
        h2.EdgeColor = 'k';
        h2.NumBins = 1;
        h2.BinWidth = 0.1;
    end
end

for i=PBS3
    if ismember(i, flagged)==0
        h3=histogram(B(i,:))
        h3.FaceColor = 'c';
        h3.EdgeColor = 'k';
        h3.NumBins = 1;
        h3.BinWidth = 0.1;
    end
end
xlabel('beta')
ylabel('Counts')
title('Distribution of beta coefficients, 2 min=red, 8 min=blue, 16 min=cyan')
saveas(gcf, 'hist_B_PBS.png')
saveas(gcf, 'hist_B_PBS.fig')

close all

%% Histogram plots of PBS C
flagged=find(A>1.5|A<0);
figure(4), hold on
for i=PBS1
    if ismember(i, flagged)==0
        h1=histogram(C(i,:))
        h1.FaceColor = 'r';
        h1.EdgeColor = 'k';
        h1.NumBins = 1;
        h1.BinWidth = 0.1;
    end
end
xlabel('C')
ylabel('Counts')
title('Distribution of C coefficients, 2 min=red')
saveas(gcf, 'hist_C_PBS1.png')
saveas(gcf, 'hist_C_PBS1.fig')

figure(5), hold on
for i=PBS2
    if ismember(i, flagged)==0
        h2=histogram(C(i,:))
        h2.FaceColor = 'b';
        h2.EdgeColor = 'k';
        h2.NumBins = 1;
        h2.BinWidth = 0.1;
    end
end
xlabel('C')
ylabel('Counts')
title('Distribution of C coefficients,8 min=blue')
saveas(gcf, 'hist_C_PBS2.png')
saveas(gcf, 'hist_C_PBS2.fig')

figure(6), hold on
for i=PBS3
    if ismember(i, flagged)==0
        h3=histogram(C(i,:))
        h3.FaceColor = 'c';
        h3.EdgeColor = 'k';
        h3.NumBins = 1;
        h3.BinWidth = 0.1;
    end
end
xlabel('C')
ylabel('Counts')
title('Distribution of C coefficients, 16 min=cyan')
saveas(gcf, 'hist_C_PBS3.png')
saveas(gcf, 'hist_C_PBS3.fig')

close all

%% Histogram plots of PBS delta
flagged=find(A>1.5|A<0);
figure(4), hold on
for i=PBS1
    if ismember(i, flagged)==0
        h1=histogram(D(i,:))
        h1.FaceColor = 'r';
        h1.EdgeColor = 'k';
        h1.NumBins = 1;
        h1.BinWidth = 0.1;
    end
end
title('Distribution of delta coefficients, 2 min=red')
saveas(gcf, 'hist_D_PBS1.png')
saveas(gcf, 'hist_D_PBS1.fig')

figure(5), hold on
for i=PBS2
    if ismember(i, flagged)==0
        h2=histogram(D(i,:))
        h2.FaceColor = 'b';
        h2.EdgeColor = 'k';
        h2.NumBins = 1;
        h2.BinWidth = 0.1;
    end
end
xlabel('delta')
ylabel('Counts')
title('Distribution of delta coefficients, 8 min=blue')
saveas(gcf, 'hist_D_PBS2.png')
saveas(gcf, 'hist_D_PBS2.fig')

figure(6), hold on
for i=PBS3
    if ismember(i, flagged)==0
        h3=histogram(D(i,:))
        h3.FaceColor = 'c';
        h3.EdgeColor = 'k';
        h3.NumBins = 1;
        h3.BinWidth = 0.1;
    end
end
xlabel('delta')
ylabel('Counts')
title('Distribution of delta coefficients, 16 min=cyan')
saveas(gcf, 'hist_D_PBS3.png')
saveas(gcf, 'hist_D_PBS3.fig')

%close all
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