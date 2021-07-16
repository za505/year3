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
basename='07142021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07142021_analysis/' basename '/' basename '_colony1/' basename '_phase2/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07142021_analysis/' basename '/' basename '_colony1/'  basename '_mNeonGreen/' basename '_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07142021_analysis/' basename '/' basename '_colony1/' basename '_mNeonGreen/' basename '_aligned']}; 
recrunch=0;
replot=0;
troubleshoot=0;
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
    icell_green=nan(ncells, T);
    norm_green=nan(ncells, T);
    f={};
    
    for i=1:length(channels)
        
        cd(channels{i});
        %intensity=nan(ncells, T);
        
        for t=1:T
             
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            for n=1:ncells
                icell_green(n,t)=mean(im(pixels{n,t}));
            end
            
        end
        
        for n=1:ncells
            norm_green(n,:) = icell_green(n,:)./icell_green(n,1);
        end    
        
        for n=1:ncells
            
            flag=0;
            guess=1;
            
            while flag==0
                
                [f{n}, answer]=tryPlot(time, norm_green, n, guess);
                guess=answer;
                
                if answer==1
                    cd(savedir)
                    savePlot(time, norm_green, n, f{n}, basename);
                    flag=1;
                end
            end
        end
            

    end
    
end


%% Plot data
if replot==1
    cd(savedir)
    %plot to see single traces of mNeonGreen cells
    figure(1), hold on
    for i=1:height(icell_green)
        plot(time, icell_green(i,:))
        x=time(end);
        y=icell_green(i,end);
        text(x,y, num2str(i));
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('Cellular Intensity of mNeonGreen vs Time')
    xlabel('Time (s)')
    ylabel('Cellular Intensity (A.U.)')
    saveas(gcf, [basename,'_fullGreen.fig'])
    saveas(gcf, [basename,'_fullGreen.png'])

    figure(2), hold on
    for n=1:height(icell_green)
        plot(time, lcell(n, :))
        x=time(end);
        y=lcell(n,end);
        text(x,y, num2str(n));
    end
    %xline(tpt, '--', {'Membrane Lysis'})
    title('mNeonGreen Cell Length vs Time')
    xlabel('Time (s)')
    ylabel('Length (\mum)')
    saveas(gcf, [basename,'_LTGreen.fig'])
    saveas(gcf, [basename,'_LTGreen.png'])

%     %plot to see single traces of mNeonGreen cells
%     for i=1:height(icell_green)
%         figure('Name', num2str(i))
%         plot(f{i}, time, icell_green(i,:))
%         title(['Cellular Intensity of mNeonGreen vs Time, ' '#' num2str(i)])
%         xlabel('Time (s)')
%         ylabel('Cellular Intensity (A.U.)')
%         saveas(gcf, [basename '_' num2str(i) '_fitGreen.fig'])
%         saveas(gcf, [basename '_' num2str(i) '_fitGreen.png'])
%         close
%     end
end

%% Troubleshooting
if troubleshoot==1
    cd(channels{1}); 
     for t=1:T
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
end
%% Calculate time constant histograms
% dir1=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07142021_analysis/' basename '_colony1/' basename '_GFP/' basename '_figures'];
% dir2=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07142021_analysis/' basename '_colony2/' basename '_GFP/' basename '_figures'];
% 
% cd(dir1)
% load([basename '_dm.mat'], 'f')
% f1=f;
% 
% cd(dir2)
% load([basename '_dm.mat'], 'f')
% f2=f;
% 
% tconst=[];
% for i=1:2
%     if i==1
%         for n=1:length(f1)
%             tconst=[tconst; -1/f1{n}.b];
%         end
%     elseif n==2
%         for n=1:length(f2)
%             tconst=[tconst; -1/f2{n}.b];
%         end
%     end
% end
% 
% cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07142021_analysis/'])
% bins=length(f1)+length(f2);
%  histogram(tconst,bins)
%  title('Histogram of Time Constants, n=37')
%  xlabel('Time Scale Tau (s-)')
%  saveas(gcf, [basename '_histogram.fig'])
%  saveas(gcf, [basename '_histogram.png'])

cd(savedir)
save([basename '_dm.mat'])
    
%% Functions
function savePlot(time, norm_green, i, newModel, basename) 
                figure('Name', num2str(i))
                plot(newModel, time, norm_green(i,:))
                title(['Cellular Intensity of mNeonGreen vs Time, ' '#' num2str(i)])
                xlabel('Time (s)')
                ylabel('Cellular Intensity (A.U.)')
                saveas(gcf, [basename '_' num2str(i) '_fitGreen.fig'])
                saveas(gcf, [basename '_' num2str(i) '_fitGreen.png'])
                close  
end

function [newModel, answer]= tryPlot(time, norm_green, n, guess)
            
            if guess==1
                
                [xData, yData] = prepareCurveData(time, norm_green(n,:));
                newModel=fit(xData, yData, 'exp1');
                plot(newModel, time, norm_green(n,:))
                pause
                
                prompt = 'Is this a good fit? 1=Yes, 2=Try Two-Term Exponential, 3=Try Linear, 4=Try Quadratic, 5=Try Piecewise Linear, 6=Piecewise Cubic ';
                answer = input(prompt)
                
                close
                
            elseif guess==2
                
                [xData, yData] = prepareCurveData(time, norm_green(n,:));
                newModel=fit(xData, yData, 'exp2');
                plot(newModel, time, norm_green(n,:))
                pause 
                
                prompt = 'Is this a good fit? 1=Yes, 2=Try Two-Term Exponential, 3=Try Linear, 4=Try Quadratic, 5=Try Piecewise Linear, 6=Piecewise Cubic ';
                answer = input(prompt)
                
                close 
                
            elseif guess==3
                
                [xData, yData] = prepareCurveData(time, norm_green(n,:));
                newModel=fit(xData, yData, 'poly1');
                plot(newModel, time, norm_green(n,:))
                pause
                
                prompt = 'Is this a good fit? 1=Yes, 2=Try Two-Term Exponential, 3=Try Linear, 4=Try Quadratic, 5=Try Piecewise Linear, 6=Piecewise Cubic ';
                answer = input(prompt)
                
                close
                
            elseif guess==4
                
                [xData, yData] = prepareCurveData(time, norm_green(n,:));
                newModel=fit(xData, yData, 'poly2');
                plot(newModel, time, norm_green(n,:))
                pause
                
                prompt = 'Is this a good fit? 1=Yes, 2=Try Two-Term Exponential, 3=Try Linear, 4=Try Quadratic, 5=Try Piecewise Linear, 6=Piecewise Cubic ';
                answer = input(prompt)
                
                close
                
            elseif guess==5
                
                [xData, yData] = prepareCurveData(time, norm_green(n,:));
                newModel=fit(xData, yData, 'linearinterp');
                plot(newModel, time, norm_green(n,:))
                pause
                
                prompt = 'Is this a good fit? 1=Yes, 2=Try Two-Term Exponential, 3=Try Linear, 4=Try Quadratic, 5=Try Piecewise Linear, 6=Piecewise Cubic ';
                answer = input(prompt)

                close
                
            elseif guess==6
                
                [xData, yData] = prepareCurveData(time, norm_green(n,:));
                newModel=fit(xData, yData, 'cubicinterp');
                plot(newModel, time, norm_green(n,:))
                pause
                
                prompt = 'Is this a good fit? 1=Yes, 2=Try Two-Term Exponential, 3=Try Linear, 4=Try Quadratic, 5=Try Piecewise Linear, 6=Piecewise Cubic ';
                answer = input(prompt)
                
                close
                
            else
                
                prompt = 'Are you sure? 0=No, 1=Yes '
                response = input(prompt)
                
                if response==1
                    [xData, yData] = prepareCurveData(time, norm_green(n,:));
                    newModel=fit(xData, yData, 'exp1');
                    plot(newModel, time, norm_green(n,:))
                    pause(0.3)

                    message='default: exponential'
                    answer=1;
                else
                    answer=2;
                end
                
            end
            
end
        