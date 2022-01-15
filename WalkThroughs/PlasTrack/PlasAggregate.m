%PlasAggregate.m
%Purpose: to aggregate and analysis plasmolysis data quantified through
%PlasTrack.m

clear, close all
%% User Input
ponA2_exp = {'12162021_Exp1', '12162021_Exp3'}; %, '12172021_Exp1', '12182021_Exp1'};
WT_exp = {'12162021_Exp2', '12162021_Exp4'}; %, '12172021_Exp2', '12182021_Exp2'};
maindir=['/Users/zarina/Documents/MATLAB/MatlabReady/Msmegmatis_analysis/01142022'];
%% Analysis: Do dPonA2 GFP cells plasmolyze more than WT GFP cells?
%pre-allocate variables
ponA2_dir = [];
WT_dir = [];

ponA2_prePerc = [];
ponA2_postPerc = [];

WT_prePerc = [];
WT_postPerc = [];

ponA2_preD={};
ponA2_postD={};

WT_preD={};
WT_postD={};

ponA2_imPre=[];
ponA2_imPost=[];

WT_imPre=[];
WT_imPost=[];

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
    
    %first, the pre and post plasm. perc. as well as the midline variables
    cd(maindir)
    load(ponA2_dir(i).name, 'pre_perc', 'post_perc', 'pre_percDist', 'post_percDist', 'imGpre', 'imGpost');
    ponA2_prePerc = [ponA2_prePerc, pre_perc'];
    ponA2_postPerc = [ponA2_postPerc, post_perc'];
    
    ponA2_preD=[ponA2_preD; pre_percDist];
    ponA2_postD=[ponA2_postD; post_percDist];
       
    ponA2_imPre(:, :, i)=imGpre;
    ponA2_imPost(:, :, i)=imGpost;
    
end

for i=1:length(WT_dir)
    
    %first, the pre and post plasm. perc. as well as the directory where
    %the mline (midline) variable is stored
    cd(maindir)
    load(WT_dir(i).name, 'pre_perc', 'post_perc', 'pre_percDist', 'post_percDist', 'imGpre', 'imGpost');
    WT_prePerc = [WT_prePerc, pre_perc'];
    WT_postPerc = [WT_postPerc, post_perc'];
    
    WT_preD=[WT_preD; pre_percDist];
    WT_postD=[WT_postD; post_percDist];
    
    WT_imPre(:, :, i)=imGpre;
    WT_imPost(:, :, i)=imGpost;
    
end

%plot the data
cd(maindir)
% labels = {'WT Pre-shock', 'ponA2 Pre-shock', 'WT Post-shock', 'ponA2 Post-shock'};
% xval = [zeros(1,length(WT_prePerc)), ones(1, length(ponA2_prePerc)), repelem(2, length(WT_postPerc)), repelem(3, length(ponA2_postPerc))];
% yval = [WT_prePerc, ponA2_prePerc, WT_postPerc, ponA2_postPerc];
% 
% figure
% scatter(xval, yval)
% title('WT vs ponA2 Plasmolysis')
% ylabel('Percent Plasmolysis')
% xticklabels(labels)
% xticks([0 1 2 3])
% xlim([-0.5 3.5])
% ylim([0 100])
% saveas(gcf, ['WTvsPonA2.png'])
% saveas(gcf, ['WTvsPonA2.fig'])
% 
% %are the pre-shock plamsolysis percentages stat. sig. different b/t WT and
% %ponA2?
% ttest2(WT_prePerc, ponA2_prePerc) %0=fail to reject the null
% %are the post-shock plamsolysis percentages stat. sig. different b/t WT and
% %ponA2?
% ttest2(WT_postPerc, ponA2_postPerc) %1=reject the null
% %are the pre-shock and post-shock plamsolysis percentages stat. sig.
% %different for WT?
% ttest2(WT_prePerc, WT_postPerc) %1=reject the null
% %are the pre-shock and post-shock plamsolysis percentages stat. sig.
% %different for ponA2?
% ttest2(ponA2_prePerc, ponA2_postPerc) %1=reject the null
% 
% perc_mat = [WT_prePerc(1:161)', ponA2_prePerc(1:161)', WT_postPerc(1:161)', ponA2_postPerc(1:161)'];
% [p, tbl, stats] = anova1(perc_mat, labels);
% saveas(gcf, ['WTvsPonA2_anova.png'])
% saveas(gcf, ['WTvsPonA2_anova.fig'])
% c = multcompare(stats)

%is there a difference in the pixel distribution of ponA2 pre and post
%shock cells? There should be
% figure, hold on
% histogram(cell2mat(ponA2_preD))
% histogram(cell2mat(ponA2_postD))
% xlim([0 50])
% xlabel('Percent Cell Length')
% ylabel('Number of Pixels (ponA2)')
% legend({'pre-shock', 'post-shock'})
% 
% figure, hold on
% histogram(cellfun(@median, ponA2_preD), height(ponA2_preD))
% histogram(cellfun(@median, ponA2_postD), height(ponA2_postD))
% xlim([0 50])
% xlabel('Median Percent Cell Length')
% ylabel('Number of Cells (ponA2)')
% legend({'pre-shock', 'post-shock'})
% 
% % is there a difference in the pixel distribution of WT pre and post
% % shock cells? There should be
% figure, hold on
% histogram(cell2mat(WT_preD))
% histogram(cell2mat(WT_postD))
% xlim([0 50])
% xlabel('Percent Cell Length')
% ylabel('Number of Pixels (WT)')
% legend({'pre-shock', 'post-shock'})
% 
% figure, hold on
% histogram(cellfun(@median, WT_preD), height(WT_preD))
% histogram(cellfun(@median, WT_postD), height(WT_postD))
% xlim([0 50])
% xlabel('Median Percent Cell Length')
% ylabel('Number of Cells (WT)')
% legend({'pre-shock', 'post-shock'})
% 
% % is there a difference in the pixel distribution of WT and ponA2? There
% % shouldn't be in the post-shock, according to the anova
% figure, hold on
% histogram(cell2mat(ponA2_preD))
% histogram(cell2mat(WT_preD))
% xlim([0 50])
% xlabel('Percent Cell Length')
% ylabel('Number of Pixels Pre-shock')
% legend({'ponA2', 'WT'})
% 
% figure, hold on
% histogram(cell2mat(ponA2_postD))
% histogram(cell2mat(WT_postD))
% xlim([0 50])
% xlabel('Percent Cell Length')
% ylabel('Number of Pixels Post-shock')
% legend({'ponA2', 'WT'})
% 
% figure, hold on
% histogram(cellfun(@median, ponA2_preD), height(ponA2_preD))
% histogram(cellfun(@median, WT_preD), height(WT_preD))
% xlim([0 50])
% xlabel('Median Percent Cell Length')
% ylabel('Number of Cells Pre-Shock')
% legend({'ponA2', 'WT'})
% 
% figure, hold on
% histogram(cellfun(@median, ponA2_postD), height(ponA2_postD))
% histogram(cellfun(@median, WT_postD), height(WT_postD))
% xlim([0 50])
% xlabel('Median Percent Cell Length')
% ylabel('Number of Cells Post-Shock')
% legend({'ponA2', 'WT'})

%conclusion: I don't see much of a difference in any case

ponA2_preR=plasRegion(ponA2_preD);
ponA2_postR=plasRegion(ponA2_postD);

WT_preR=plasRegion(WT_preD);
WT_postR=plasRegion(WT_postD);

save(['msmeg_analysis'])

% %% Functions
% function [plasm_region]=plasRegion(percDist)
%     idx=cellfun(@isempty, percDist);
%     idx=setdiff(1:height(percDist), idx);
%     
%     percDist=percDist(idx, 1);
%     plasm_region=cell(size(percDist));
%     
%     for n=1:height(percDist)
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
%     
%     %y = discretize(x,[0 .25 .75
%     %1],'categorical',{'small','medium','large'}); maybe use this next
%     %time?
%     plasm_region=cellfun(@categorical, plasm_region, 'UniformOutput', false);
% end