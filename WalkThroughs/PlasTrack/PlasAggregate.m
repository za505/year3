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

WT_plasRegion=[];
WT_plasTotal=[];
WT_cellTotal=[];

ponA2_plasRegion=[];
ponA2_plasTotal=[];
ponA2_cellTotal=[];

%load data
for i=1:length(ponA2_dir)
    
    %load the pre and post count plasmolysis data
    %note: 1=polar, 2=subpolar, 3=midcell
    cd(maindir)
    load(ponA2_dir(i).name, 'preCount', 'postCount');
    ponA2_preCount = [ponA2_preCount; preCount];
    ponA2_postCount = [ponA2_postCount; postCount];
    [ponA2_plasRegion(i, :), ponA2_plasTotal(i, :), ponA2_cellTotal(i, :)] = plasRegion(ponA2_postCount);

end

for i=1:length(WT_dir)
    
    cd(maindir)
    load(WT_dir(i).name, 'preCount', 'postCount');
    WT_preCount = [WT_preCount; preCount];
    WT_postCount = [WT_postCount; postCount];
    [WT_plasRegion(i, :), WT_plasTotal(i, :), WT_cellTotal(i, :)] = plasRegion(WT_postCount);
end

%assign the data into groups based on regions where plasmolysis was
%observed, row 1 = polar, row 2 = subpolar, row 3 = midcell
% [ponA2_plasRegion, ponA2_plasTotal, ponA2_cellTotal] = plasRegion(ponA2_postCount);
% [WT_plasRegion, WT_plasTotal, WT_cellTotal] = plasRegion(WT_postCount);

%plot the data

%Plasmolysis Probability/Cell
cd(maindir)
figure, hold on
b=bar([mean(ponA2_plasRegion./ponA2_cellTotal,1)', mean(WT_plasRegion./WT_cellTotal,1)'])
err1=std(ponA2_plasRegion./ponA2_cellTotal,0, 1);
err2=std(WT_plasRegion./WT_cellTotal,0, 1);
errorbar(b(1).XEndPoints, b(1).YData, err1, 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k')
errorbar(b(2).XEndPoints, b(2).YData, err2, 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k')
legend({'ponA2', 'WT'})
ylabel('Plasmolysis Probability/Cell', 'FontSize', 15)
xticklabels({'Polar', '', 'Subpolar', '', 'Midcell'})
xlabel('Region', 'FontSize', 15)

%Average Number of Plasmolysis Bays/Cell
cd(maindir)
figure, hold on
b1=bar(1, mean(ponA2_plasTotal./ponA2_cellTotal,1))
b2=bar(2, mean(WT_plasTotal./WT_cellTotal,1))

err1=std(ponA2_plasTotal./ponA2_cellTotal,0, 1);
err2=std(WT_plasTotal./WT_cellTotal,0, 1);
errorbar(b1.XEndPoints, b1.YData, err1, 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k')
errorbar(b2.XEndPoints, b2.YData, err2, 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k')

legend({'ponA2', 'WT'})
ylabel('Number of Plasmolysis Bays/Cell', 'FontSize', 15)
xticklabels({''})
xlabel('Strain', 'FontSize', 15)
% saveas(gcf, 'plasProbCell.png')
% saveas(gcf, 'plasProbCell.fig')

A=ponA2_plasTotal./ponA2_cellTotal;
B=WT_plasTotal./WT_cellTotal;
C=[A,B];

b(1).YData=A;
b(2).YData=B;

err1=std(ponA2_plasRegion./ponA2_cellTotal,0, 1);
err2=std(WT_plasRegion./WT_cellTotal,0, 1);
errorbar(b(1).XEndPoints, b(1).YData, err1, 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k')
errorbar(b(2).XEndPoints, b(2).YData, err2, 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k')

ylabel('# plasmolysis bays/cell', 'FontSize', 15)
xlabel('')
%legend({'ponA2', 'WT'})
xticklabels({'', 'ponA2', '', 'WT'})
% saveas(gcf, 'plasPer.png')
% saveas(gcf, 'plasPer.fig')

A=ponA2_plasTotal;
B=WT_plasTotal;
C=[A,B];

b(1,1).YData=A;
b(1,2).YData=B;
ylabel('Total Incidence of Plasmolysis', 'FontSize', 15)
xlabel('')
%legend({'ponA2', 'WT'})
xticklabels({'', 'ponA2', '', 'WT'})
% saveas(gcf, 'plasTotal.png')
% saveas(gcf, 'plasTotal.fig')
% 
% save(['01272022_msmegAnalysis'])

%% Functions
function [plasm_region, plasm_total, cell_total]=plasRegion(plasCount)
    idx=cellfun(@isempty, plasCount);
    idx=setdiff(1:height(plasCount), idx);
    
    plasCount=plasCount(idx, 1);
    plasm_region=zeros(3,1);
    cell_total=height(plasCount);
    
    for n=1:height(plasCount)
        if length(intersect(plasCount{n, 1}, 1:3))==1
            if plasCount{n,1}==1
                plasm_region(1,1)=plasm_region(1,1)+1;
            elseif plasCount{n,1}==2
                plasm_region(2,1)=plasm_region(2,1)+1;
            elseif plasCount{n,1}==3
                plasm_region(3,1)=plasm_region(3,1)+1;
            end
        elseif length(intersect(plasCount{n, 1}, 1:3))>1
            if sum(intersect(plasCount{n, 1}, 1:3))==3
                plasm_region(1,1)=plasm_region(1,1)+1;
                plasm_region(2,1)=plasm_region(2,1)+1;
            elseif sum(intersect(plasCount{n, 1}, 1:3))==4
                plasm_region(1,1)=plasm_region(1,1)+1;
                plasm_region(3,1)=plasm_region(3,1)+1;
            elseif sum(intersect(plasCount{n, 1}, 1:3))==5
                plasm_region(2,1)=plasm_region(2,1)+1;
                plasm_region(3,1)=plasm_region(3,1)+1;
            elseif sum(intersect(plasCount{n, 1}, 1:3))==5
                plasm_region(1,1)=plasm_region(1,1)+1;                
                plasm_region(2,1)=plasm_region(2,1)+1;
                plasm_region(3,1)=plasm_region(3,1)+1;
            end
        else
            continue   
        end
    end   
    
    %plasm_region=plasm_region';
    plasm_total=sum(plasm_region);
    
    %y = discretize(x,[0 .25 .75
    %1],'categorical',{'small','medium','large'}); maybe use this next
    %time?
%     plasm_region=cellfun(@categorical, plasm_region, 'UniformOutput', false);
end