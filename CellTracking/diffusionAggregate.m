%diffusionAggregate.m
%Author: Zarina Akbary
%Date: 04/15/2022
%Purpose: this is basically a sandbox for me to look through my data and
%keep track of how I do it. 

clear, close all

%% the goal is to aggregate the 20% intensity data
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/raw/';
dirsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/aggregate/';
basename = '05062022_Exp1';

cd(dirpath)
    
datadir = dir([basename '*']);
[intensity, bgintensity, adjintensity, lcell, time] = dataAggregate(datadir);

cd(dirsave)
save(basename)

%% Functions
%to aggregate data from a given experiment
function [intensity, bgintensity, adjintensity, lCell, time]=dataAggregate(datadir)
    
    %pre-allocate variables
    intensity=[];
    bgintensity=[];
    adjintensity=[];
    lCell=[];
        
    %go through the data for each position
    for i=1:length(datadir)

        %load decayMeasure .mat file
        cd(datadir(i).folder)
        load(datadir(i).name, '-regexp', 'intensity$', 'time', 'lcell')
        
        %check to see if the background exists
        if exist('bg_intensity', 'var')
            bgintensity=[bgintensity; bg_intensity];
            adj_intensity=icell_intensity-bg_intensity;
        else
            bgintensity=NaN;
            adjintensity=NaN;
        end
        
        %concatenate intensity and length data
        intensity = [intensity; icell_intensity];
        adjintensity = [adjintensity; adj_intensity];
        lCell = [lCell; lcell];
        
    end
         
end