%BTColorDecode.m
%3/28/19, Rico Rojas
%
%To be used after BacTrack.  
%Decodes multiple strains based on fluorescent labels.
%Save color decode images as an image sequence in a separate directory with
%phase image first.

clear, close all

%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basename='05262021_Exp1'; legend_str={'GFP','TADA'}; thresholds=[600 20000];  %t1=[15 27 39 51 63];  t2=[17 29 41 53 66];%Plasmolysis lysis of MG1655, KB242, EHEC pooled, occlusion microscopy
%basename='032819_5', legend_str={'EHEC','MG1655','KB242'},thresholds=[2000 1200], t1=[13 28 40 52 64];  t2=[15 30 42 55 67];;%Plasmolysis lysis of MG1655, KB242, EHEC pooled, occlusion microscopy
%basename='032819_7',legend_str={'MG1655','KB242','EHEC'},thresholds=[400 1000], t1=[13 25 37 49 60];  t2=[15 27 39 52 63];%Plasmolysis lysis of MG1655, KB242, EHEC pooled, occlusion microscopy

ncolors=3;%Number of channels including dark cells

movie_directory=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05262021_analysis/' basename '/' basename '_phase/' basename '_erased'];%Directory in which movie is saved
color_directory=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05262021_analysis/' basename];%Directory in which color 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


curdir=cd;
cd(color_directory); 
fluo_directory=dir('*.tif');
cd(movie_directory); 
phase_directory=dir('*.tif');

cd([color_directory '/' basename '_phase/' basename '_figures']); 
load([basename '_BTphase'],'T','pixels','v','lcell','tmid','ncells','B')

T=2;
pixels=pixels(:,1:T);
B=B(:,1:T);
v=v(:,1:T-1);
lcell=lcell(:, 1:T);
tmid=tmid(:,1:T);


    cd(color_directory); 
    imagename=fluo_directory(1).name;
    TestIm=imread(imagename);
    cd(movie_directory); 
    imagename=phase_directory(T).name;
    RefIm=imread(imagename);
    
    [output NewImFT]=dftregistration(fft2(RefIm),fft2(TestIm),10);
    NewIm=abs(ifft2(NewImFT));
    shft(1)=output(3);
    shft(2)=output(4);

    vav_col=zeros(ncolors,T-1);
    
    
    cd(color_directory); 
for i=1:ncolors-1
    
    imagename=fluo_directory(i).name;
    im=imread(imagename);
    
    im=imtranslate(im,shft);
    
    intensities_temp=zeros(ncells,1);
        for j=1:ncells
            intensities_temp(j)=mean(im(pixels{j,T}));
        end
    icell{i}=intensities_temp;
    
    positive_cells{i}=find(intensities_temp>thresholds(i));    
    
    v_col{i}=v(positive_cells{i},:);
    vav_col(i,:)=nanmean(v_col{i});
   
    figure,hist(nonzeros(intensities_temp)),pause
   
    figure,
    imshow(im,[])
    %imo=imoverlay(im,ed10,[1 0 0]);
    %imshow(imo) 
    hold on
   for k=1:length(positive_cells{i})
    %   plot(P{k}(:,2),P{k}(:,1),'-y')
      plot(B{positive_cells{i}(k),T}(:,1),B{positive_cells{i}(k),T}(:,2),'-r')
   end
  pause
end

final_cells=find(~isnan(lcell(:,T)));
colored_cells=union(positive_cells{2},positive_cells{3});
positive_cells{1}=setdiff(final_cells,colored_cells);

v_col{1}=v(positive_cells{1},:);
vav_col(1,:)=nanmean(v_col{1});

Leff_col=zeros(size(vav_col));
for i=1:ncolors
    Leff_col(i,:)=EffectiveLength(tmid,vav_col(i,:));
end

figure,hold on, title(basename)
for i=1:ncolors
    plot(tmid,Leff_col(i,:))
end
xlabel('Time (s)')
ylabel('Effective Population Averaged Length (/mum)')
legend(legend_str)
fig2pretty
% 
% Strn_col=(Leff_col(:,t1)-Leff_col(:,t2))./Leff_col(:,t2);
% 
% shock_mag=[0 50 100 200 400];
% shock_mag_int=[0:1:400];
% cmap=colormap(jet);
% figure,hold on,title(basename)
% for i=1:ncolors
%     X = [ones(length(shock_mag),1) shock_mag'];
%     b = X\Strn_col(i,:)';
%     strn_int=shock_mag_int*b(2)+b(1);
%     stiffness(i)=1./b(2);
%     p(2*i-1)=plot(shock_mag, Strn_col(i,:)*100,'ok','MarkerSize',12,'MarkerFaceColor',cmap(round(i*64/ncolors),:));
%     p(2*i)=plot(shock_mag_int,strn_int*100,':','Color',cmap(round(i*64/ncolors),:));
%     set(get(get(p(2*i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% end
% legend(p(1:2:end),legend_str)
% ylabel('Strain (%)')
% xlabel('Shock Magnitude (mM)')
% fig2pretty
% 
% figure,hold on
% title(basename)
% bar(stiffness)
% ylabel('Envelope Stiffness')
% xticks([1:3])
% xticklabels(legend_str)
% fig2pretty

cd(curdir)

save([basename '_BTColDec'])