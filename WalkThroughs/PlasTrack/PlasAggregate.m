%PlasAggregate.m
%Purpose: to aggregate and analysis plasmolysis data quantified through
%PlasTrack.m

clear, close all
%% User Input
ponA2_exp = {'12162021_Exp1', '12162021_Exp3'}; %, '12172021_Exp1', '12182021_Exp1'};
WT_exp = {'12162021_Exp2', '12162021_Exp4'}; %, '12172021_Exp2', '12182021_Exp2'};
maindir=['/Users/zarina/Documents/MATLAB/MatlabReady/Msmegmatis_analysis/01142022'];
%% Analysis: Do dPonA2 GFP cells plasmolyze more/in different regions than WT GFP cells?
%pre-allocate variables
ponA2_dir = [];
WT_dir = [];

ponA2_preCount = {};
ponA2_postCount = {};

WT_preCount = {};
WT_postCount = {};

ponA2_postA=[];
WT_postA=[];

%set up structure with file names
cd(maindir)
for i=1:length(ponA2_exp)
    ponA2_dir = [ponA2_dir; dir([ponA2_exp{i} '*'])];
end

for i=1:length(WT_exp)
    WT_dir = [WT_dir; dir([WT_exp{i} '*'])];
end

%load data
for i=1:length(ponA2_dir)
    
    %load the pre and post count plasmolysis data
    %note: 1=polar, 2=subpolar, 3=midcell
    cd(maindir)
    load(ponA2_dir(i).name, 'preCount', 'postCount');
    ponA2_preCount = [ponA2_preCount; preCount];
    ponA2_postCount = [ponA2_postCount; postCount];
    
end

for i=1:length(WT_dir)
    
    cd(maindir)
    load(WT_dir(i).name, 'preCount', 'postCount');
    WT_preCount = [WT_preCount; preCount];
    WT_postCount = [WT_postCount; postCount];
    
end

%assign the data into groups based on regions where plasmolysis was
%observed
%note: the groups are as follows, 1=no plasmolysis bays, 2=only polar,
%3=polar and subpolar, 4=polar and midcell, 5=subpolar, 6=subpolar and
%midcell, 7=midcell, 8=polar, subpolar, and midcell
ponA2_postA = plasRegion(ponA2_postCount);
WT_postA = plasRegion(WT_postCount);

%tabluate the data to get a table of frequencies
ponA2_postTab=tabulate(ponA2_postA);
WT_postTab=tabulate(WT_postA);

%both matrices should have the same number of rows
if height(ponA2_postTab)<8 & height(WT_postTab)<8
    mv=setdiff(1:8, 1:height(ponA2_postTab));
    amend = [mv', zeros(length(mv), 1), zeros(length(mv), 1)];
    ponA2_postTab=[ponA2_postTab; amend];
    
    mv=setdiff(1:8, 1:height(WT_postTab));
    amend = [mv', zeros(length(mv), 1), zeros(length(mv), 1)];
    WT_postTab=[WT_postTab; amend];
elseif height(ponA2_postTab)<8
    mv=setdiff(1:8, 1:height(ponA2_postTab));
    amend = [mv', zeros(length(mv), 1), zeros(length(mv), 1)];
    ponA2_postTab=[ponA2_postTab; amend];
elseif height(WT_postTab)<8
    mv=setdiff(1:8, 1:height(WT_postTab));
    amend = [mv', zeros(length(mv), 1), zeros(length(mv), 1)];
    WT_postTab=[WT_postTab; amend];
end

%plot the data
cd(maindir)
figure, hold on
bar([ponA2_postTab(:, 2), WT_postTab(:, 2)])
legend({'ponA2', 'WT'})
ylabel('Counts', 'FontSize', 15)
xticklabels('')
saveas(gcf, 'plasDistribution.png')
saveas(gcf, 'plasDistribution.fig')

save(['msmeg_analysis'])

%% Functions
function [plasm_region]=plasRegion(plasCount)
    idx=cellfun(@isempty, plasCount);
    idx=setdiff(1:height(plasCount), idx);
    
    plasCount=plasCount(idx, 1);
    plasm_region=nan(height(plasCount),1);
    
    for n=1:height(plasCount)
        if length(intersect(plasCount{n, 1}, 1:3))==1
            if plasCount{n,1}==1
                plasm_region(n,1)=2;
            elseif plasCount{n,1}==2
                plasm_region(n,1)=5;
            elseif plasCount{n,1}==3
                plasm_region(n,1)=7;
            end
        elseif length(intersect(plasCount{n, 1}, 1:3))>1
            if sum(intersect(plasCount{n, 1}, 1:3))==3
                plasm_region(n,1)=3;
            elseif sum(intersect(plasCount{n, 1}, 1:3))==4
                plasm_region(n,1)=4;
            elseif sum(intersect(plasCount{n, 1}, 1:3))==5
                plasm_region(n,1)=6;
            elseif sum(intersect(plasCount{n, 1}, 1:3))==5
                plasm_region(n,1)=8;
            end
        else
            plasm_region(n,1)=1;    
        end
    end
    
%         for p=1:height(percDist{n,1})
%             pdist=percDist{n,1}(p,1);
%             if pdist <12.5
%                plasm_region{n,1}(p,1)="polar";
%             elseif pdist >=12.5 & pdist <25
%                 plasm_region{n,1}(p,1)="subpolar";
%             elseif pdist >=25
%                 plasm_region{n,1}(p,1)="mid-cell";
%             end
%         end
%     end
    
    %y = discretize(x,[0 .25 .75
    %1],'categorical',{'small','medium','large'}); maybe use this next
    %time?
%     plasm_region=cellfun(@categorical, plasm_region, 'UniformOutput', false);
end