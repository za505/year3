%decayMeasure.m
%Zarina Akbary, updated 02/09/22
%Calculates changes in fluor. intensity. Incorporates BTfluo.m code.

clear, close all

%INSTRUCTIONS FOR USE:
%run BacTrack.m first

%INPUT
%basename: experiments of interest
%dirname: where .mat files are stored
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%icell_intensity = cell x time matrix of fluoresc. intensities
%bg_intensity = 1 x time matrix of mean background intensities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='10232021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2022_Spring/10232021_analysis/' basename '/' basename '_colony2/' basename '_phase/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2022_Spring/02102022_reanalysis'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2022_Spring/10232021_analysis/' basename '/' basename '_colony2/' basename '_mNeonGreen/' basename '_aligned']}; 
recrunch=0;
replot=1;
troubleshoot=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savedir)
    load([basename '_colony2_dm.mat'])
    replot=0;
    troubleshoot=2;
else
    
    for i=1:length(channels)
        cd(channels{i}); 
        fluo_directory{i}=dir('*.tif');
    end
    
    %go to directory where .mat files are stored
    cd(dirname)
    load([basename '_BTphase'], 'B', 'T', 'ncells', 'time', 'pixels', 'lcell')

    %pre-allocate variables
    icell_intensity=nan(ncells, T);
    bg_intensity=nan(1,T);
    septa=zeros(ncells, 1);
    time=time./60;
    
    %find the final pre-lysis frame 
    dl=diff(lcell, 1, 2);
    lvg=mean(dl, 1, 'omitnan');
    [~, lidx]=min(lvg);
    lidx=lidx+2;   
    
    for i=1:length(channels)
        
        cd(channels{i});
        %midx=round(T/2);
        
        imagename=fluo_directory{i}(lidx).name;
        [p1, p2]=getBackground(imagename);
        
        for t=1:T
            t
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            bg_intensity(1, t)=measureBackground(im, p1, p2);

            for n=1:ncells
                icell_intensity(n,t)=mean(im(pixels{n,t}));
            end
            
        end  
       
    end    
end
%% Troubleshooting
if troubleshoot==1
    cd(channels{1}); 
     for t=1:T
            t
            imagename=fluo_directory{1}(1).name;
            im=imread(imagename);


           figure
           imshow(im, [])
           hold on
           for k=1:ncells
                if isempty(B{k,t})==0
                   plot(B{k,t}(:,1),B{k,t}(:,2),'-g')
                else
                    continue
                end
           end
          pause
          close all
     end
     
elseif troubleshoot==2
    cd(channels{1});
     for k=1:ncells
            k
            imagename=fluo_directory{1}(lidx).name;
            im=imread(imagename);


           figure
           imshow(im, [])
           hold on
           for t=1:T
                if isempty(B{k,t})==0
                    plot(B{k,t}(:,1),B{k,t}(:,2),'-g')
                else
                    continue
                end
           end
          pause
          close all
     end
end

prompt = 'Which cells have septa? ';
idx = input(prompt)
septa(idx, 1) = 1; 

%% Plot data
if replot==1
    
    cd(savedir)
    
    %plot cellular fluorescence traces
        figure(1), hold on
    for i=1:height(icell_intensity)
        if septa(i,1)
            plot(time, icell_intensity(i,:), '-m')
        else
            plot(time, icell_intensity(i,:), '-g')
        end
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Intensity vs Time')
    %subtitle('blue = halved','Color','blue')
    xlabel('Time (min)')
    ylabel('Cellular Intensity (A.U.)')
    saveas(gcf, [basename,'_fullIntensity_dm.fig'])
    saveas(gcf, [basename,'_fullIntensity_dm.png'])

    %plot background fluorescence traces
    figure(2), hold on
    plot(time, bg_intensity, '-b')
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Intensity vs Time')
    %subtitle('blue = halved','Color','blue')
    xlabel('Time (min)')
    ylabel('Background Intensity (A.U.)')
    saveas(gcf, [basename,'_bgIntensity_dm.fig'])
    saveas(gcf, [basename,'_bgIntensity_dm.png'])
    
%     figure(2), hold on
%     for i=1:height(adj_intensity)
%         plot(time, adj_intensity(i,:), '-g')
%     end
%     %xline(tpt, '--', {'Membrane Lysis'})
%     title('Adjusted Intensity of phase vs Time')
%     %subtitle('blue = halved','Color','blue')
%     xlabel('Time (min)')
%     ylabel('Cellular Intensity (A.U.)')
%     saveas(gcf, [basename,'_adjIntensity_dm.fig'])
%     saveas(gcf, [basename,'_adjIntensity_dm.png'])
%     
%     %plot to see single traces of phase cells
%     figure(3), hold on
%     for i=1:height(norm_intensity)
%         plot(time, norm_intensity(i,:), '-g')
%     end
%     %xline(tpt, '--', {'Membrane Lysis'})
%     title('Normalized Intensity of phase vs Time')
%     %subtitle('blue = halved','Color','blue')
%     xlabel('Time (min)')
%     ylabel('Cellular Intensity (A.U.)')
%     saveas(gcf, [basename,'_normIntensity_dm.fig'])
%     saveas(gcf, [basename,'_normIntensity_dm.png'])
% 
%     figure(4), hold on
%     for n=1:height(norm_intensity)
%         plot(time, lcell(n, 1:end), '-g')
%     end
%     %xline(tpt, '--', {'Membrane Lysis'})
%     title('Cell Length vs Time')
%     %subtitle('blue = halved','Color','blue')
%     xlabel('Time (min)')
%     ylabel('Length (\mum)')
%     saveas(gcf, [basename,'_LTintensity_dm.fig'])
%     saveas(gcf, [basename,'_LTintensity_dm.png'])
% 
end

cd(savedir)
save([basename '_colony2_dm.mat'])

% cd('/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/02072022_analysis/MatFiles/')
% save([basename '_colony2_dm.mat'])

%% Functions
 function [p1, p2]=getBackground(imagename)

         %Load last image
         %imagename=fluo_directory{i}(t).name;
         im2=imread(imagename);

         %Determine Background
         figure,imshow(im2,[]), hold on, title('Select Background')
         k=waitforbuttonpress;
         set(gcf,'Pointer')
         hold on
         axis manual
         point1=get(gca,'CurrentPoint');
         finalRect=rbbox;
         point2=get(gca,'CurrentPoint');
         point1=point1(1,1:2);
         point2=point2(1,1:2);
         point1(point1<1)=1;
         point2(point2<1)=1;
         p1=min(point1,point2);%Calculate locations
         p2=max(point1,point2);
         offset = abs(point1-point2);%And dimensions
         x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
         y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
         plot(x,y)
         p1=round(p1);
         p2=round(p2);  
 end 

 function bglevel = measureBackground(im2, p1, p2)

         %Load last image
         %imagename=fluo_directory{i}(t).name;
         %im2=imread(imagename);

         %Determine background
         backim=im2(p1(2):p2(2),p1(1):p2(1));
         [counts,bins]=imhist(backim);
         [~,binnum]=max(counts);
         maxpos=bins(binnum);
         bglevel=mean(mean(backim));

 end 