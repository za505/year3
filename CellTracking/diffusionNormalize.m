%Author: Zarina Akbary
%Date: 03/26/2022
%Purpose: To normalize diffusion data from dm.mat scripts

clear, close all

%directory guide
%untreated, 1 minute frame rate, 100% intensity = '10232021_Exp1' & '10262021_Exp1' 
%untreated, 1 minute frame rate, 20% intensity = '02212022_Exp2' 
%untreated, 5 minute frame rate, 100% intensity = '02122022_Exp1' 
%untreated, 10 minute frame rate, 100% intensity = '02122022_Exp2'
%untreated, 20 minute frame rate, 100% intensity = '02092022_Exp1'
%untreated, 20 minute frame rate, 20% intensity = '02212022_Exp1'
 
%2 minute PBS incubation, 1 minute frame rate, 100% intensity = '11192021_Exp2' & '12082021_Exp3' 
%20 minute PBS incubation, 1 minute frame rate, 100% intensity = '11192021_Exp1' & '11302021_Exp1'
%60 minute PBS incubation, 1 minute frame rate, 100% intensity = '10232021_Exp2' & '10262021_Exp2'
%120 minute PBS incubation, 1 minute frame rate, 100% intensity = '01142022_Exp1'

%LB + 20 mM Mg2+, 1 minute frame rate, 100% intensity = '01172022_Exp1' 
%LB + 10 mM EDTA, 1 minute frame rate, 100% intensity = '01172022_Exp2' 
%LB + 0.5 ug/mL tunicamycin, 1 minute frame rate, 100% intensity = '01242022_Exp1' 
%LB + 1 ug/mL vancomycin, 1 minute frame rate, 100% intensity = '01262022_Exp1'
%spent LB, 10 minute frame rate, 100% intensity = '02122022_Exp3' & '02192022_Exp1'

%untreated, 1.2 s frame rate, 100% intensity = '11202021_Exp1' 
%untreated, 2 s frame rate, 100% intensity = '12082021_Exp1'
%untreated, 3 s frame rate, 100% intensity = '12082021_Exp2'

%untreated, 1.76 s frame rate, 20% intensity = '02192022_Exp2' 
%untreated, 2.3 s frame rate, 20% intensity = '02192022_Exp3'
%untreated, 3 s frame rate, 20% intensity = '02192022_Exp4'

%% User Input

%location and names of of the *.mat files 
dirsave='/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis';
basenames={'10232021_Exp1', '10262021_Exp1', '02212022_Exp2', '02122022_Exp1','02122022_Exp2', '02092022_Exp1', '02212022_Exp1', '11192021_Exp2', '12082021_Exp3', '11192021_Exp1', '11302021_Exp1', '10232021_Exp2', '10262021_Exp2', '01142022_Exp1', '01172022_Exp1', '01172022_Exp2', '01242022_Exp1', '01262022_Exp1', '02122022_Exp3', '02192022_Exp1','11202021_Exp1', '12082021_Exp1', '12082021_Exp2', '02192022_Exp2', '02192022_Exp3', '02192022_Exp4'};

%% normalize data
imstarts=[6, 5, 6, 6, 6, 6, 6, 10, 10, 4, 4, 3, 5, 4, 6, 6, 3, 5, 6, 6, NaN, NaN, NaN, NaN, NaN, NaN]; %based on time vector
idxes={[2:4], [2:6], [2:4], [2:3], [2:3], [2:3], [2:3], [2:4], [2:5], [2:4], [2:4], [2:4], [2:4], [], [2:7], [2:4], [2:4], [], [2:3], [2:3], [], [], [], [], [], []}; %based on tme vector

for b=1:length(basenames)

    basename=basenames{b};
    intp=idxes{b};
    
    cd([dirsave '/rawFiles'])
    
    datadir=dir([basename '*']);
    imstart=imstarts(b);
    
    if ismember(b, [21:26]) %controls
        
        intensity=[];
        bgintensity=[];
        adjintensity=[];
        normintensity=[];
        lcell=[];
        
        for j=1:height(datadir)
            cd([dirsave '/rawFiles']);
            load(datadir(j).name);
            
            
            [intensity1, bgintensity1, adjintensity1, normintensity1, lCell, time, tme, imstart]=controlsNormalize(datadir);
            intensity=[intensity;intensity1];
            bgintensity=[bgintensity; bgintensity1];
            adjintensity=[adjintensity; adjintensity1];
            normintensity=[normintensity; normintensity1];
            lcell=[lcell; lCell];
            
        end
        
        cd([dirsave '/normalizedFiles'])
        save([basename '_norm.mat'], 'intensity', 'bgintensity', 'adjintensity', 'normintensity', 'lcell', 'imstart', 'time', 'tme');
    
    else
        
        intensity=[];
        bgintensity=[];
        adjintensity=[];
        normintensity=[];
        lcell=[];
        
        for j=1:height(datadir)
            cd([dirsave '/rawFiles'])
            load(datadir(j).name);
            
            
            [intensity1, bgintensity1, adjintensity1, normintensity1, lCell, time, tme, imstart]=dataNormalize(datadir, imstart, intp);
            intensity=[intensity;intensity1];
            bgintensity=[bgintensity; bgintensity1];
            adjintensity=[adjintensity; adjintensity1];
            normintensity=[normintensity; normintensity1];
            lcell=[lcell; lCell];
            
        end
        
        cd([dirsave '/normalizedFiles'])
        save([basename '_norm.mat'], 'intensity', 'adjintensity', 'normintensity', 'lcell', 'imstart', 'time', 'tme', 'intp');
    end
    
end

%% Functions
%to aggregate data and normalize
function [intensity, bgintensity, adjintensity, normintensity, lCell, time, tme, imstart]=dataNormalize(datadir, imstart, idx)
    
    %pre-allocate variables
    intensity=[];
    bgintensity=[];
    lCell=[];
        
    %go through the data for each position
    for i=1:length(datadir)

        %load decayMeasure .mat file
        cd(datadir(i).folder)
        load(datadir(i).name, '-regexp', 'intensity$', 'time', 'lcell')
        
        
        %we only need the variables icell_intensity, time, lcell, and--if
        %it exists--bg_intensity
        
        if exist('bg_intensity', 'var')
            bgintensity=[bgintensity; bg_intensity];
        else
            bgintensity=NaN;
        end
        
        [nrow, ~]=size(icell_intensity);
        
        for n=1:nrow
            if ~isnan(icell_intensity(n, imstart))
                intensity=[intensity; icell_intensity(n, :)];
                lCell=[lCell; lcell(n, :)];
            end
        end

        if i==1
            tme=time; %pre-set new time vector
        end

    end 
       
    %subtract final fluor. value from the trace
    adjintensity=intensity-intensity(:, end); 
    adjintensity(adjintensity<0)=NaN;      
        
    %adjust the time points in the fluor matrix and initialize to first
    %frame
    normintensity=adjintensity(:, imstart:end)./adjintensity(:, imstart);
    [nrow, ncol]=size(normintensity);
    
    %interpolate the fluor values during detergent perfusion 
    for n=1:height(normintensity)        
            if ~isempty(idx)
                x=setdiff(1:ncol, idx); %x=time
                v=normintensity(n, x); %v=intensity
                vq=interp1(x, v, idx); %vq=interpolation at query time pts
                normintensity(n,idx)=vq;
            end
    end
    
    %adjust the time vector
    tme=tme(imstart:end)-tme(imstart);     
        
end

%to aggregate controls data and normalize
function [intensity, bgintensity, adjintensity, normintensity, lCell, time, tme, imstart]=controlsNormalize(datadir)

    %pre-allocate variables
    intensity=[];
    lCell=[];
    bgintensity=[];
        
    %go through the data for each position
    for i=1:length(datadir)

        %load decayMeasure .mat file
        cd(datadir(i).folder)
        load(datadir(i).name, 'icell_intensity', 'bg_intensity', 'time', 'lcell')
        
        intensity=[intensity; icell_intensity];
        lCell=[lCell; lcell];
        bgintensity=[bgintensity; bg_intensity];

        if i==1
            tme=time; %pre-set new time vector
        end

    end
        

    %find the initial post-lysis frame by identifying the first
    %time point where dt matches the final dt
    dt=round(diff(time));
    if dt(1)==dt(end)
%         dl=diff(lCell, 1);
%         lvg=mean(dl, 1, 'omitnan');
%         imstart=find(lvg<0, 1, 'first')+2;
        imstart=1;
        
        %remove cells without a value for imstart
        idx=find(~isnan(intensity(:, imstart)));
        intensity=intensity(idx,:);
        lCell=lCell(idx,:);    

        %subtract final fluor. value from the trace and set limit of
        %detection
        adjintensity=intensity-bgintensity; 
        adjintensity(adjintensity<200)=NaN; 
    
    else
        if dt(end)<1 %discrepancies in dt might be found in the meta data, this is the ad hoc fix
            imstart=find(dt==dt(end), 1, 'first')+2;
        else 
            imstart=find(dt==dt(end), 1, 'first');
        end
    
  
    %remove cells without a value for imstart
    idx=find(~isnan(intensity(:, imstart)));
    intensity=intensity(idx,:);
    lCell=lCell(idx,:);    
        
    %subtract final fluor. value from the trace and set limit of
    %detection
    adjintensity=intensity-intensity(:, end); 
    adjintensity(adjintensity<200)=NaN; 
    
    end
        
    %find when all values are below the limit of detection
    [~, ncol]=size(adjintensity);
    nsum=sum(isnan(adjintensity));
    if max(nsum)==ncol
        imend=min(find(nsum==ncol));
    else 
        imend=ncol;
    end
    
    %adjust the time points in the fluor matrix and initialize to first
    %frame
    normintensity=adjintensity(:, imstart:imend)./adjintensity(:, imstart);

    %adjust the time vector
    tme=tme(imstart:imend)-tme(imstart);     
        
end
