%sanity_checks.m
%Author: Zarina Akbary
%Date: 04/15/2022
%Purpose: this is basically a sandbox for me to look through my data and
%keep track of how I do it. 

clear, close all

%Is there a difference in the background/autofluorescence value of cells
%perfused with LB vs LB + LB + IPTG? Let's look at cases where the frame
%rate is 1 minute from beginning to end.

%% initial perfusion with LB + IPTG, final perfusion with LB
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/rawFiles/';
basename = '04042022_Exp1';

cd(dirpath)
    
datadir = dir([basename '*']);
[intensity, bgintensity, lcell, time] = dataAggregate(datadir);

%in a plot of background intensity vs time, it doesn't appear as though the
%rate of photobleaching changes. But this is just by eye. There must be a
%more rigorous approach.
figure
plot(time, bgintensity)

%ignore time points when detergent is perfused
time1 = time;
time1(5:9) = NaN;

bgTrace1 = bgintensity;
bgTrace1(:, 5:9) = NaN;

%% initial perfusion with LB + IPTG, final perfusion with LB + IPTG. Includes the CY5 channel.
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/rawFiles/';
basename = '10232021_Exp1';

cd(dirpath)
    
datadir = dir([basename '*']);
[intensity, bgintensity, lcell, time] = dataAggregate(datadir);

%in a plot of background intensity vs time, it doesn't appear as though the
%rate of photobleaching changes. But this is just by eye. There must be a
%more rigorous approach.
figure
plot(time, bgintensity)

%ignore time points when detergent is perfused
time2 = time;
time2(5:9) = NaN;

bgTrace2 = bgintensity;
bgTrace2(:, 5:9) = NaN;

%% initial perfusion with LB + IPTG, final perfusion with LB + IPTG. Includes the CY5 channel.
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/rawFiles/';
basename = '10262021_Exp1'

cd(dirpath)
    
datadir = dir([basename '*']);
[intensity, bgintensity, lcell, time] = dataAggregate(datadir);

%in a plot of background intensity vs time, it doesn't appear as though the
%rate of photobleaching changes. But this is just by eye. There must be a
%more rigorous approach.
figure
plot(time, bgintensity)

%ignore time points when detergent is perfused
time3 = time;
time3(5:9) = NaN;

bgTrace3 = bgintensity;
bgTrace3(:, 5:9) = NaN;

%% initial perfusion with LB + IPTG, final perfusion with LB. 20% intensity
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/rawFiles/';
basename = '02212022_Exp2';

cd(dirpath)
    
datadir = dir([basename '*']);
[intensity, bgintensity, lcell, time] = dataAggregate(datadir);

%in a plot of background intensity vs time, it doesn't appear as though the
%rate of photobleaching changes. But this is just by eye. There must be a
%more rigorous approach.
figure
plot(time, bgintensity)

%ignore time points when detergent is perfused
time4 = time;
time4(5:9) = NaN;

bgTrace4 = bgintensity;
bgTrace4(:, 5:9) = NaN;

%% plot to compare the three
figure, hold on
plot(time1, bgTrace1, '-r')
plot(time2, bgTrace2, '-b')
plot(time3, bgTrace3, '-g')
plot(time4, bgTrace4, '-k')

%% Functions
%to aggregate data from a given experiment
function [intensity, bgintensity, lCell, time]=dataAggregate(datadir)
    
    %pre-allocate variables
    intensity=[];
    bgintensity=[];
    lCell=[];
        
    %go through the data for each position
    for i=1:length(datadir)

        %load decayMeasure .mat file
        cd(datadir(i).folder)
        load(datadir(i).name, '-regexp', 'intensity$', 'time', 'lcell')

        %concatenate intensity and length data
        intensity = [intensity; icell_intensity];
        lCell = [lCell; lcell];
        
        %check to see if the background exists
        if exist('bg_intensity', 'var')
            bgintensity=[bgintensity; bg_intensity];
        else
            bgintensity=NaN;
        end
        
    end
         
end