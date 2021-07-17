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
basename='07162021_Exp4';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07162021_analysis/' basename '/' basename '_colony4/' basename '_phase/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07162021_analysis/' basename '/' basename '_colony4/'  basename '_mNeonGreen/' basename '_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07162021_analysis/' basename '/' basename '_colony4/' basename '_mNeonGreen/' basename '_aligned']}; 
recrunch=0;
replot=1;
troubleshoot=2;
tidx=9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savedir)
    load([basename '_dm.mat'])
else
    
    for i=1:length(channels)
        cd(channels{i}); 
        fluo_directory{i}=dir('*.tif');
    end
    
    %go to directory where .mat files are stored
    cd(dirname)
    load([basename '_BTphase'], 'B', 'T', 'ncells', 'time', 'pixels', 'lcell')

    %pre-allocate variables
    icell_green=nan(ncells, T-tidx+1);
    norm_green=nan(ncells, T-tidx+1);
    time=time(tidx:end)./60;
    
    for i=1:length(channels)
        
        cd(channels{i});
        %intensity=nan(ncells, T);
        
        for t=tidx:T
            j=t-(tidx-1);
            
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            for n=1:ncells
                icell_green(n,j)=mean(im(pixels{n,t}));
            end
            
        end
        
        %normalize the intensity traces
        for n=1:ncells
            norm_green(n,:) = icell_green(n,:)./icell_green(n,1);
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
                    plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
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
                    plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
                else
                    continue
                end
           end
          pause
          close all
     end
end

%% Does only half the cell lose fluor?
%halfie=nan(ncells, 1);
n_halfie=[1,2,3,5,6,9,12,16,17];
halfie=ismember([1:ncells], n_halfie)';
% cd(channels{1}); 
% 
% for k=1:ncells
%     
%     %tL=max(find(cellfun(@isempty, B(k,:))==0)); %where is the last time point with cells?
%     
% %     v = VideoWriter(strcat('mNeonGreen_', num2str(k), '_intensity'),'MPEG-4');
% %     open(v);
%     
%     for t=T
%         
%         imagename=fluo_directory{1}(t).name;
%         im=imread(imagename);
% 
%         figure
%         imshow(im, [])
%         hold on
%        
%         if isempty(B{k,t})==0
%             plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
%         end
%         
% %         frame = getframe(gcf);
% %         writeVideo(v,frame);
% %         pause(1)
% %         clf
% 
%         pause
%     end
%      
% %      close(v)
% %      close all
%      
%      prompt = 'Is only half the cell fluor or dark? 0=No, 1=Yes ';
%      answer = input(prompt)
%      
%      halfie(k)=answer;
%      
%  end

%% combine these variables into a table
dataTable=table(lcell, icell_green, norm_green, halfie, 'VariableNames', {'cell length', 'intensity', 'normalized intensity','halfie'});

%% Plot data
if replot==1
    
    cd(savedir)
    
    %plot to see single traces of mNeonGreen cells
    figure(1), hold on
    for i=1:height(norm_green)
        if halfie(i)==0
            plot(time, norm_green(i,:), '-r')
            x=time(end);
            y=norm_green(i,end);
            text(x,y, num2str(i));
        elseif halfie(i)==1
            plot(time, norm_green(i,:), '-b')
            x=time(end);
            y=norm_green(i,end);
            text(x,y, num2str(i));
        end
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Normalized Intensity of mNeonGreen vs Time')
    subtitle('blue = halved','Color','blue')
    xlabel('Time (min)')
    ylabel('Cellular Intensity (A.U.)')
    saveas(gcf, [basename,'_normGreen.fig'])
    saveas(gcf, [basename,'_normGreen.png'])

    figure(2), hold on
    for n=1:height(norm_green)
        if halfie(n)==0
            plot(time, lcell(n, tidx:end), '-r')
            x=time(end);
            y=lcell(n,end);
            text(x,y, num2str(n));
        elseif halfie(n)==1
            plot(time, lcell(n, tidx:end), '-b')
            x=time(end);
            y=lcell(n,end);
            text(x,y, num2str(n));
        end
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Cell Length vs Time')
    subtitle('blue = halved','Color','blue')
    xlabel('Time (min)')
    ylabel('Length (\mum)')
    saveas(gcf, [basename,'_LTGreen.fig'])
    saveas(gcf, [basename,'_LTGreen.png'])

    figure(3), hold on
    for i=1:height(icell_green)
        if halfie(i)==0
            plot(time, icell_green(i,:), '-r')
            x=time(end);
            y=icell_green(i,end);
            text(x,y, num2str(i));
        elseif halfie(i)==1
            plot(time, icell_green(i,:), '-b')
            x=time(end);
            y=icell_green(i,end);
            text(x,y, num2str(i));
        end
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Intensity of mNeonGreen vs Time')
    subtitle('blue = halved','Color','blue')
    xlabel('Time (min)')
    ylabel('Cellular Intensity (A.U.)')
    saveas(gcf, [basename,'_fullGreen.fig'])
    saveas(gcf, [basename,'_fullGreen.png'])
    
%     %plot to see single traces of mNeonGreen cells
%     for i=1:height(icell_green)
%         figure('Name', num2str(i))
%         plot(f{i}, time, icell_green(i,:))
%         title(['Cellular Intensity of mNeonGreen vs Time, ' '#' num2str(i)])
%         xlabel('Time (min)')
%         ylabel('Cellular Intensity (A.U.)')
%         saveas(gcf, [basename '_' num2str(i) '_fitGreen.fig'])
%         saveas(gcf, [basename '_' num2str(i) '_fitGreen.png'])
%         close
%     end
end

cd(savedir)
save([basename '_dm.mat'])
    