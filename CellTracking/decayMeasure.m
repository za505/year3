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
basename='10232021_Exp2';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/10232021_analysis/' basename '/' basename '_p002_colony6/' basename '_p002_phase/' basename '_p002_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/10232021_analysis/'  basename '/' basename '_p002_colony6/'  basename '_p002_mNeonGreen/' basename '_p002_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Fall/10232021_analysis/' basename '/' basename '_p002_colony6/' basename '_p002_mNeonGreen/' basename '_p002_aligned']}; 
recrunch=0;
replot=1;
troubleshoot=0;
tidx=8; %the first fluor image or the first post-lysis image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savedir)
    load([basename '_colony6_dm3.mat'])
else
    
    for i=1:length(channels)
        cd(channels{i}); 
        fluo_directory{i}=dir('*.tif');
    end
    
    %go to directory where .mat files are stored
    cd(dirname)
    load([basename '_BTphase'], 'B', 'T', 'ncells', 'time', 'pixels', 'lcell')

    %pre-allocate variables
    if tidx==1
        icell_intensity=nan(ncells, T);
        norm_intensity=nan(ncells, T);
        time=time(tidx:end)./60;
        time=time-time(1);
    else
        icell_intensity=nan(ncells, T-tidx+1);
        norm_intensity=nan(ncells, T-tidx+1);
        time=time(tidx:end)./60;
        time=time-time(1);
    end
    
    for i=1:length(channels)
        
        cd(channels{i});
        %intensity=nan(ncells, T);
        
        for t=tidx:T
            j=t-(tidx-1);
            
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            for n=1:ncells
                icell_intensity(n,j)=mean(im(pixels{n,t}));
            end
            
        end
        
        %subtract the background from each trace
        for n=1:ncells
            adj_intensity(n,:) = icell_intensity(n,:)-icell_intensity(n,end);
        end 
        
        %normalize the intensity traces
        for n=1:ncells
            norm_intensity(n,:) = adj_intensity(n,:)./adj_intensity(n,1);
        end    
       
    end    
end
%% Troubleshooting
if troubleshoot==1
    cd(channels{1}); 
     for t=tidx:T
            t
            imagename=fluo_directory{1}(t).name;
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
           for t=tidx:T
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
    
    %plot to see single traces of mNeonGreen cells
    figure(1), hold on
    for i=1:height(norm_intensity)
        plot(time, norm_intensity(i,:), '-g')
        x=time(end);
        y=norm_intensity(i,end);
        text(x,y, num2str(i));
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Normalized Intensity of mNeonGreen vs Time')
    %subtitle('blue = halved','Color','blue')
    xlabel('Time (min)')
    ylabel('Cellular Intensity (A.U.)')
    saveas(gcf, [basename,'_normGreen_dm3.fig'])
    saveas(gcf, [basename,'_normGreen_dm3.png'])

    figure(2), hold on
    for n=1:height(norm_intensity)
        plot(time, lcell(n, tidx:end), '-g')
        x=time(end);
        y=lcell(n,end);
        text(x,y, num2str(n));
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Cell Length vs Time')
    %subtitle('blue = halved','Color','blue')
    xlabel('Time (min)')
    ylabel('Length (\mum)')
    saveas(gcf, [basename,'_LTGreen_dm3.fig'])
    saveas(gcf, [basename,'_LTGreen_dm3.png'])

    figure(3), hold on
    for i=1:height(icell_intensity)
        plot(time, icell_intensity(i,:), '-g')
        x=time(end);
        y=icell_intensity(i,end);
        text(x,y, num2str(i));
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Intensity of mNeonGreen vs Time')
    %subtitle('blue = halved','Color','blue')
    xlabel('Time (min)')
    ylabel('Cellular Intensity (A.U.)')
    saveas(gcf, [basename,'_fullGreen_dm3.fig'])
    saveas(gcf, [basename,'_fullGreen_dm3.png'])
    
    figure(4), hold on
    for i=1:height(adj_intensity)
        plot(time, adj_intensity(i,:), '-g')
        x=time(end);
        y=adj_intensity(i,end);
        text(x,y, num2str(i));
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Adjusted Intensity of mNeonGreen vs Time')
    %subtitle('blue = halved','Color','blue')
    xlabel('Time (min)')
    ylabel('Cellular Intensity (A.U.)')
    saveas(gcf, [basename,'_adjGreen_dm3.fig'])
    saveas(gcf, [basename,'_adjGreen_dm3.png'])
    

end

cd(savedir)
save([basename '_colony6_dm3.mat'])

cd('/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/10272021_analysis')
save([basename '_colony6_dm3.mat'])