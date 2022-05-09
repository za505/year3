%matlab_gui.m
%Purpose: to make a gui

clear, close all

%% re-analysis of mNeonGreen data
%I want to see if I can improve the fit for all the traces
%data=icell_green;

%% user input
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Downloads/NYU/Year3_2021_Fall/presentations/09062021";

data = struct('experiment', [], 'colony', struct('norm_green', [], 'time', [], 'a', [], 'b', [], 'c', [], 'd', [], 'fitModel', []));
cd(datadir)

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

%%
buttonPlot(data(2).colony(1).norm_green, data(2).colony(1).time)

function buttonPlot(data, time)

% Create a figure window
fig= uifigure;

%Give the figure a 3-by-2 grid
g = uigridlayout(fig);
g.RowHeight = {35, 35, 35, 35, 35, 35, 35,35, '1x'}; %the height of the top two rows is 35 px while the bottom row height is variable
g.ColumnWidth = {150,'1x'}; %the width of the first column is 150 px while the second column is variable

%create drop-down menu to choose which cell intensity trace to fit the data
%to
items=[1:height(data)];
%pickLabel = uilabel(g,'Text','Cell:');
pickCell = uidropdown(g, 'Items', string(items), 'ItemsData', items, 'ValueChangedFcn', @(pickCell,event) selection(pickCell,data));
%pickCell.ItemsData=items;
pickCell.Layout.Row = 1;
pickCell.Layout.Column = 2;

% Create a UI axes
ax = uiaxes(g); 
ax.Layout.Row = [2 9];
ax.Layout.Column = 2;

% Create a buttons to fit data
btn1 = uibutton(g,'push', 'Text', 'Single-Term Exponential', ...
               'ButtonPushedFcn', @(btn1,event) plotExp1(btn1,ax, data(pickCell.Value,:), time));
btn1.Layout.Row = 1;
btn1.Layout.Column = 1;

btn2 = uibutton(g,'push', 'Text', 'Two-Term Exponential', ...
               'ButtonPushedFcn', @(btn2,event) plotExp2(btn2,ax, data(pickCell.Value,:), time));
btn2.Layout.Row = 2;
btn2.Layout.Column = 1;

btn3 = uibutton(g,'push', 'Text', 'Linear', ...
               'ButtonPushedFcn', @(btn3,event) plotLinear(btn3,ax, data(pickCell.Value,:), time));
btn3.Layout.Row = 3;
btn3.Layout.Column = 1;

btn4 = uibutton(g,'push', 'Text', 'Quadratic', ...
               'ButtonPushedFcn', @(btn4,event) plotQuadratic(btn4,ax, data(pickCell.Value,:), time));
btn4.Layout.Row = 4;
btn4.Layout.Column = 1;

btn5 = uibutton(g,'push', 'Text', 'Piecewise Linear', ...
               'ButtonPushedFcn', @(btn5,event) plotPLinear(btn5,ax, data(pickCell.Value,:), time));
btn5.Layout.Row = 5;
btn5.Layout.Column = 1;

btn6 = uibutton(g,'push', 'Text', 'Piecewise Quadratic', ...
               'ButtonPushedFcn', @(btn6,event) plotPQuadratic(btn6,ax, data(pickCell.Value,:), time));
btn6.Layout.Row = 6;
btn6.Layout.Column = 1;

btn7 = uibutton(g,'push', 'Text', 'Piecewise Cubic', ...
               'ButtonPushedFcn', @(btn7,event) plotPCubic(btn7,ax, data(pickCell.Value,:), time));
btn7.Layout.Row = 7;
btn7.Layout.Column = 1;

btn8 = uibutton(g,'push', 'Text', 'Save', ...
               'ButtonPushedFcn', @(btn8,event) savePlot(btn8,ax, data(pickCell.Value,:), time));
btn8.Layout.Row = 8;
btn8.Layout.Column = 1;
        
end

% Create ValueChangedFcn callback:
function [cellData]=selection(pickCell, data)
    val=pickCell.Value;
    cellData = data(val, :);
end

% Create the function for the ButtonPushedFcn callback
function plotExp1(btn1,ax, cellData, time)
        [xData, yData] = prepareCurveData(time, cellData);
        p=fit(xData, yData, 'exp1')
        scatter(ax, xData,yData), hold on
        plot(p, xData, yData)
end

function plotExp2(btn2,ax)
        x = linspace(0,2*pi,100);
        y = sin(x);
        plot(ax,x,y)
end

function plotLinear(btn3,ax)
        x = linspace(0,2*pi,100);
        y = sin(x);
        plot(ax,x,y)
end

function plotQuadratic(btn4,ax)
        x = linspace(0,2*pi,100);
        y = sin(x);
        plot(ax,x,y)
end

function plotPLinear(btn5,ax)
        x = linspace(0,2*pi,100);
        y = sin(x);
        plot(ax,x,y)
end

function plotPQuadratic(btn6,ax)
        x = linspace(0,2*pi,100);
        y = sin(x);
        plot(ax,x,y)
end

function plotPCubic(btn7,ax)
        x = linspace(0,2*pi,100);
        y = sin(x);
        plot(ax,x,y)
end

function savePlot(btn8,ax)
        x = linspace(0,2*pi,100);
        y = sin(x);
        plot(ax,x,y)
end

%% first, I want to plot how intensity changes overtime
% %the change when switching from PBS + 5% NLS to *PBS [+ Mg2]
% colony1=['/Users/zarina/Documents/MATLAB/MatlabReady/06062021_analysis/06062021_Exp1_colony1'];
% colony2=['/Users/zarina/Documents/MATLAB/MatlabReady/06062021_analysis/06062021_Exp1_colony2'];
% 
% cd(colony1)
% load('06062021_Exp1_dm.mat', 'icell_green', 'icell_mcherry', 'nGreen', 'nCherry')
% colony1_green=icell_green;
% colony1_mcherry=icell_mcherry;
% 
% cd(colony2)
% load('06062021_Exp1_dm.mat', 'icell_green', 'icell_mcherry', 'nGreen', 'nCherry', 'time', 'tpt', 'tidx')
% colony2_green=icell_green;
% colony2_mcherry=icell_mcherry;
% 
% combo_green=[colony1_green; colony2_green];
% combo_mcherry=[colony1_mcherry; colony2_mcherry];
% 
% n_green=height(combo_green);
% n_mcherry=height(combo_mcherry);
% 
% expGreen=cell(n_green,1);
% expCherry=cell(n_mcherry,1);
% 
% time=time(tidx);
% cutoff=find(time==tpt);
% time2=time(cutoff:end)-tpt;
% 
% norm_green=[];
% norm_mcherry=[];
% 
% cd('..')
% 
% for i=1:n_green
%     i
%     norm_green=[norm_green;combo_green(i,cutoff:end)/combo_green(i,cutoff)];
%     [xData, yData]=prepareCurveData(time2, norm_green(i,:));
%     expGreen{i,1}=fit(xData, yData, 'exp1');
%     plot(expGreen{i,1}, xData, yData)
%     title(['Fitting mNeonGreen Fluorescence Intensity, n= ' num2str(i)]) 
%     xlabel('Time (s)')
%     ylabel('Normalized Intensity (A.U.)')
% 
%     saveas(gcf, ['06062021_Exp1_mNeonGreenFit_' num2str(i) '.png'])
%     saveas(gcf, ['06062021_Exp1_mNeonGreenFit_' num2str(i) '.fig'])
%     
%     close
% end
% 
% for i=1:n_mcherry
%     i
%     norm_mcherry=[norm_mcherry; combo_mcherry(i,cutoff:end)./combo_mcherry(i,cutoff)];
%     [xData, yData]=prepareCurveData(time2, norm_mcherry(i,:));
%     expCherry{i,1}=fit(xData, yData, 'exp1');
%     plot(expCherry{i,1}, xData, yData)
%     title(['Fitting mCherry Fluorescence Intensity, n= ' num2str(i)]) 
%     xlabel('Time (s)')
%     ylabel('Normalized Intensity (A.U.)')
%     
%     saveas(gcf, ['06062021_Exp1_mCherryFit_' num2str(i) '.png'])
%     saveas(gcf, ['06062021_Exp1_mCherryFit_' num2str(i) '.fig'])
%     
%     close
% end
% 
% poorfit_green=[2,7,10,16,18,19];
% poorfit_mcherry=[3,7,9,10,11,12];
% 
% tconst_green=[];
% tconst_mcherry=[];
% 
% figure,hold on
% for i=1:n_green
%     if ismember(i, poorfit_green)==0
%         plot(time2./60, norm_green(i,:))
%         %tconst_green=[tconst_green -1/expGreen{i}.b];
%     end
% end
% title('mNeonGreen Fluorescence Intensity Traces') 
% xlabel('Time (min)')
% ylabel('Normalized Intensity (A.U.)')  
% saveas(gcf, ['06062021_Exp1_mNeonGreenTraces.png'])
% saveas(gcf, ['06062021_Exp1_mNeonGreenTraces.fig'])
% pause, close
% 
% figure,hold on
% for i=1:n_mcherry
%     if ismember(i, poorfit_mcherry)==0
%         plot(time2./60, norm_mcherry(i,:)) 
%         %tconst_mcherry=[tconst_mcherry -1/expCherry{i}.b];
%     end
% end
% title('mCherry Fluorescence Intensity Traces') 
% xlabel('Time (min)')
% ylabel('Normalized Intensity (A.U.)')  
% saveas(gcf, ['06062021_Exp1_mCherryTraces.png'])
% saveas(gcf, ['06062021_Exp1_mCherryTraces.fig'])
% pause, close
% 
% histogram(tconst_mcherry./60, length(tconst_mcherry))
% title('time constant (tau) for mCherry')
% xlabel('tau (min-)')
% saveas(gcf, ['06062021_Exp1_mCherryHistogram.png'])
% saveas(gcf, ['06062021_Exp1_mCherryHistogram.fig'])
% pause, close
% 
% histogram(tconst_green./60, length(tconst_green))
% title('time constant (tau) for mNeonGreen')
% xlabel('tau (min-)')
% saveas(gcf, ['06062021_Exp1_greenHistogram.png'])
% saveas(gcf, ['06062021_Exp1_greenHistogram.fig'])
% pause, close