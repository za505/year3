%decayMeasure.m
%Zarina Akbary, updated 07/04/21
%Calculates changes in fluor. intensity. Incorporates BTfluo.m code.

clear, close all

%INSTRUCTIONS FOR USE:
%run BacTrack.m first

%INPUT
%basename: experiments of interest
%dirname: where .mat files are stored
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%timescale=array of time scales
%Ts=array of T (number of frames)


%switchF=array of switchFrames
%A=array of A values (A*e^(-alpha*t)+y0)
%f=cell of coeff for exponential eqxn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='02052022_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2022_Spring/02052022_analysis/'  basename '/' basename '_colony5/' basename '_phase/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2022_Spring/02052022_analysis/'  basename '/' basename '_colony5/' basename '_mNeonGreen/' basename '_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2022_Spring/02052022_analysis/'  basename '/' basename '_colony5/' basename '_mNeonGreen/' basename '_aligned']}; 
recrunch=0;
replot=1;
troubleshoot=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savedir)
    load([basename '_colony5_dm.mat'])
    replot=0;
    troubleshoot=2;
else
    
    for i=1:length(channels)
        cd(channels{i}); 
        fluo_directory{i}=dir('*.tif');
    end
    
    %go to directory where .mat files are stored
    cd(dirname)
    load([basename '_BT'], 'B', 'T', 'ncells', 'time', 'pixels', 'lcell')

    %pre-allocate variables
    icell_intensity=nan(ncells, T);
%     adj_intensity=nan(ncells, T);
%     norm_intensity=nan(ncells, T);
    time=time./60;
    
    for i=1:length(channels)
        
        cd(channels{i});
        
        for t=1:T
            t
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            for n=1:ncells
                icell_intensity(n,t)=mean(im(pixels{n,t}));
            end
            
        end
        
%         %subtract the background from each trace
%         for n=1:ncells
%             adj_intensity(n,:) = icell_intensity(n,:)-icell_intensity(n,end);
%         end 
%         
%         %normalize the intensity traces
%         for n=1:ncells
%             norm_intensity(n,:) = adj_intensity(n,:)./adj_intensity(n,1);
%         end    
       
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
            imagename=fluo_directory{1}(T).name;
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

%% Plot data
if replot==1
    
    cd(savedir)
    
    %plot phase fluorescence traces
        figure(1), hold on
    for i=1:height(icell_intensity)
        plot(time, icell_intensity(i,:), '-g')
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Intensity of phase vs Time')
    %subtitle('blue = halved','Color','blue')
    xlabel('Time (min)')
    ylabel('Cellular Intensity (A.U.)')
    saveas(gcf, [basename,'_fullIntensity_dm.fig'])
    saveas(gcf, [basename,'_fullIntensity_dm.png'])
    
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
save([basename '_colony5_dm.mat'])

cd('/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/01252022_analysis/MatFiles/')
save([basename '_colony5_dm.mat'])