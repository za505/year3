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
basename='06062021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_reanalysis/'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_aligned']; ['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_aligned']}; 
recrunch=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(dirname)
    load([basename '_dm.mat'])
else
    
    for i=1:length(channels)
        cd(channels{i}); 
        fluo_directory{i}=dir('*.tif');
    end
    
    %go to directory where .mat files are stored
    cd(dirname)
    load([basename '_BTphase'], 'B', 'T', 'ncells', 'time', 'pixels')

    %pre-allocate variables
    icell_green=[];
    icell_mcherry=[];
    nGreen=[];
    nCherry=[];
    
    %calculate time index
    tidx=1:2:T-1; %remember, fluor images were taken every 2 frames
    tpt=480; %time at which *PBS+5% NLS is perfused
    
    for i=1:length(channels)
        
        cd(channels{i});
        intensity=nan(ncells, length(tidx));
        
        for j=1:length(tidx)
            t=tidx(j);
             
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            for n=1:ncells
                intensity(n,j)=mean(im(pixels{n,t}));
            end
            
        end
        
        for n=1:ncells
            if i==1 & mean(intensity(n, 1:5)>20000)
                icell_green=[icell_green; intensity(n,:)];
                nGreen=[nGreen n];
            elseif i==2 & mean(intensity(n, 1:5)>20000)
                icell_mcherry=[icell_mcherry; intensity(n,:)];
                nCherry=[nCherry n];
            else
                continue
            end
        end
    end
    
    cd(dirname)
    save([basename '_dm.mat'])
end

% %% Plot data
% %plot to see single traces of mNeonGreen cells
% figure(1), hold on
% for i=1:height(icell_green)
%     plot(time(tidx), icell_green(i,:))
% end
% xline(tpt, '--', {'Membrane Lysis'})
% title('Cellular Intensity of mNeonGreen vs Time')
% xlabel('Time (s)')
% ylabel('Cellular Intensity (A.U.)')
% saveas(gcf, [basename,'_fullGreen.fig'])
% saveas(gcf, [basename,'_fullGreen.png'])
%  
% %plot to see single traces of mCherry cells
% figure(2), hold on
% for i=1:height(icell_mcherry)
%    plot(time(tidx), icell_mcherry(i,:))
% end
% xline(tpt, '--', {'Membrane Lysis'})
% title('Cellular Intensity of mCherry vs Time')
% xlabel('Time (s)')
% ylabel('Cellular Intensity (A.U.)')
% saveas(gcf, [basename,'_fullCherry.fig'])
% saveas(gcf, [basename,'_fullCherry.png'])
% 
% %% make movies
% cg=0; %green count
% cc=0; %cherry count
% 
% imgreen=cell(length(nGreen),length(tidx));
% imcherry=cell(length(nCherry),length(tidx));
% 
% for n=1:ncells
%     n
%     
%     if ismember(n, nGreen)
%         cg=cg+1;
%         
%         [~, locB]=ismember(n, nGreen);
%         cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
%         
%         m=min(find(cellfun('length', B(n,:))>0)); %find the initial boundaries
%         dx=max(B{n,m}(:,1))-min(B{n,m}(:,1)); %columns are the x direction
%         dy=max(B{n,m}(:,2))-min(B{n,m}(:,2)); %rows are the y direction
%         r1=round(min(B{n,m}(:,1))-dx/2);
%         r2=round(min(B{n,m}(:,2))-dy/2);
%         d1=max([dx dy]);
%         
%         cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
%         v = VideoWriter(strcat('mNeonGreen_', num2str(cg), '_intensity'),'MPEG-4');
%         open(v);
%         
%         for j=1:length(tidx)
%             
%             t=tidx(j)
%             cd(channels{1});
%             imagename=fluo_directory{1}(t).name;
%             im=imread(imagename);
%             
%             
%             imgreen{n,j} = imcrop(im,[r1 r2 d1*2 d1*2]);
%             imshow(imgreen{n,j})
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%             pause(1)
%             clf
% 
%         end
%         
%         close(v)
%         close all
%         
% %         [imM,imN]=size(imgreen{n,m});
% %         y=1:imM; %rows are the y direction
% %         x=1:imN; %columns are the x direction
% %         [X,Y] = meshgrid(x,y);
% %         
% %         cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
% %         v = VideoWriter(strcat('mNeonGreen_', num2str(cg), '_intensity'),'MPEG-4');
% %         open(v);
% %         
% %         for j=1:length(tidx)
% %             figure, hold on
% %             t=surf(X,Y,imgreen{n,j})
% %             %rotate3d on;
% %             t.EdgeColor = 'interp';
% %             t.FaceColor = 'interp';
% %             view(0,95)
% %             frame = getframe(gcf);
% %             writeVideo(v,frame);
% %             pause(0.1)
% %             clf
% %         end
% %         
% %         close(v)
% %         close all
% 
%     elseif ismember(n, nCherry)
%         cc=cc+1;
%         
%         [~, locB]=ismember(n, nCherry);
%         cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_figures'])
%         
%         m=min(find(cellfun('length', B(n,:))>0)); %find the initial boundaries
%         dx=max(B{n,m}(:,1))-min(B{n,m}(:,1)); %columns are the x direction
%         dy=max(B{n,m}(:,2))-min(B{n,m}(:,2)); %rows are the y direction
%         r1=round(min(B{n,m}(:,1))-dx/2);
%         r2=round(min(B{n,m}(:,2))-dy/2);
%         d1=max([dx dy]);
%         
%         cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_figures'])
%         v = VideoWriter(strcat('mCherry_', num2str(cc), '_intensity'),'MPEG-4');
%         open(v);
%         for j=1:length(tidx)
%             
%             t=tidx(j)
%             cd(channels{2});
%             imagename=fluo_directory{2}(t).name;
%             im=imread(imagename);
%             
%             imcherry{n,j} = imcrop(im,[r1 r2 d1*2 d1*2]);
%             imshow(imcherry{n,j})
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%             pause(1)
%             clf
% 
%         end
%         
%         close(v)
%         close all
%         
% %         [imM,imN]=size(imcherry{n,m});
% %         y=1:imM; %rows are the y direction
% %         x=1:imN; %columns are the x direction
% %         [X,Y] = meshgrid(x,y);
% %         
% %         cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_figures'])
% %         v = VideoWriter(strcat('mCherry_', num2str(cg), '_intensity'),'MPEG-4');
% %         open(v);
% %         
% %         for j=1:length(tidx)
% %             figure, hold on
% %             t=surf(X,Y,imcherry{n,j})
% %             %rotate3d on;
% %             t.EdgeColor = 'interp';
% %             t.FaceColor = 'interp';
% %             view(0,95)
% %             frame = getframe(gcf);
% %             writeVideo(v,frame);
% %             pause(0.1)
% %             clf
% %         end
% %         
% %         close(v)
% %         close all
%     else
%         continue
%     end
% end
%         
% cd(savedir)
% save([basename 'processed_dm.mat'])
% 
% 
% %% view tracked cells
% cd(channels{1}); 
% imagename=fluo_directory{1}(1).name;
% im=imread(imagename);
% 
% for n=6 %1:length(nCherry)
%     figure, imshow(im), hold on
%     k=nGreen(n);
% for j=1:length(tidx)
%     t=tidx(j);
%     if isempty(B{k,t})==0
%         plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
%     else
%         continue
%     end
% end
% pause
% close
% end
% %% normalize cells by their initial intensity value upon lysitims
% cutoff=find(time(tidx)==tpt);
% icell_green2=icell_green(:, cutoff:end);
% icell_mcherry2=icell_mcherry(:, cutoff:end);
% 
% cutoff=find(time==tpt);
% tidx2=tidx(tidx>=cutoff);
% time2=time(tidx2)-tpt;
% 
% skip_green=[1];
% igreen=nan(height(icell_green),1);
% figure, hold on
% for n=1:height(icell_green)
%     if ismember(n,skip_green)
%         continue
%     else
%     igreen(n)=icell_green2(n,1);
%     plot(time2, icell_green2(n,:)/igreen(n,1))
%     end
% end
% 
% skip_mcherry=[1];
% imcherry=nan(height(icell_mcherry),1);
% figure, hold on
% for n=1:height(icell_mcherry)
%     if ismember(n,skip_mcherry)
%         continue
%     else
%     imcherry(n)=icell_mcherry2(n,1);
%     plot(time2, icell_mcherry2(n,:)/imcherry(n,1))
%     end
% end

%% Extract time scales from normalized data
skip_green=[1, 15];
skip_mcherry=[2, 16];

halfie_green=[9 10 11 16 17];
halfie_mcherry=[3 4 5 12 14 15];

%normalize cells by their initial intensity value upon lysis
cutoff=find(time(tidx)==tpt);
icell_green2=icell_green(:, cutoff:end);
icell_mcherry2=icell_mcherry(:, cutoff:end);

cutoff=find(time==tpt);
tidx2=tidx(tidx>=cutoff);
time2=time(tidx2)-tpt;

icelln_green=[];
icelln_mcherry=[];

for c=1:length(channels)
    for n=1:height(icell_green)
        n
        if c==1 & ismember(double(n), skip_green)==0 & ismember(double(n), halfie_green)==0
            fidx=find(isnan(icell_green2(n,:))==0);
            if isempty(fidx)==0
                icelln_green=[icelln_green; icell_green2(n,:)/icell_green2(n,fidx(1))];
            end
        elseif c==2  & ismember(double(n), skip_mcherry)==0 & ismember(double(n), halfie_mcherry)==0
            fidx=find(isnan(icell_mcherry2(n,:))==0);
            if isempty(fidx)==0
                icelln_mcherry=[icelln_mcherry; icell_mcherry2(n,:)/icell_mcherry2(n,fidx(1))];
            end
        end
    end
end 
 
fitresult=cell(height(icelln_green),1);
gof=cell(height(icelln_green),1);

for n=1:height(icelln_green)
    
    [xData, yData] = prepareCurveData(time2, icelln_green(n,:));
    
    % Set up fittype and options.
    ft = fittype( 'exp1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [1.01151960663985 -0.000417556351326276];
    
    % Fit model to data.
    [fitresult{n,:}, gof{n,:}] = fit( xData, yData, ft, opts );
    
    % Plot fit with data.
    figure( 'Name', num2str(n));
    h = plot( fitresult{n,:}, xData, yData );
    legend( h, 'mNeonGreen Fluor Intesity (A.U.) vs. Time (s)', ['cell: ' num2str(n) ', tau: ' num2str(-1/fitresult{n,:}.b)], 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'time2', 'Interpreter', 'none' );
    ylabel( 'icelln_green', 'Interpreter', 'none' );
    grid on
    pause
    
    cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
    saveas(gcf, [basename '_' num2str(n) '_expfit.fig'])
    saveas(gcf, [basename '_' num2str(n) '_expfit.png'])
    close
    
end

% for n=1:height(icelln_green)
%     
%     if ismember(n, skip_green)==0
%     [xData, yData] = prepareCurveData(time2, icelln_green(n,:));
%     
%     % Set up fittype and options.
%     ft = fittype( 'poly1' );
% 
%     % Fit model to data.
%     [ftrlt, gf] = fit( xData, yData, ft);
%     
%     % Plot fit with data.
%     figure( 'Name', num2str(n));
%     h = plot( ftrlt, xData, yData );
%     legend( h, 'mNeonGreen Fluor Intesity (A.U.) vs. Time (s)', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
%     % Label axes
%     xlabel( 'time2', 'Interpreter', 'none' );
%     ylabel( 'icelln_green', 'Interpreter', 'none' );
%     grid on
%     pause
%     
% %     cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
% %     saveas(gcf, [basename '_' num2str(n) '_expfit.fig'])
% %     saveas(gcf, [basename '_' num2str(n) '_expfit.png'])
%     close
%     end
% end


%here, we can fit each cell intensity to a nonlinear exponential eqxn
% for c=1:length(channels)
%     for i=1:height(icelln_green)
%         if c==1 & ismember(i, skip_green)==0
%             %coeff_green(i,:)=nlinfit(time2, icelln_green(i,:), @exponential, coeff0);
%             coeff_green(i,:)=lsqcurvefit(@exponential, coeff0, time2, icelln_green(i,:));
%             yhat_green(i,:)=exponential(coeff_green(i,:), time2);
%         elseif c==2 & ismember(i, skip_mcherry)==0
%             %coeff_mcherry(i,:)=nlinfit(time2, icelln_mcherry(i,:), @exponential, coeff0);
%             coeff_mcherry(i,:)=lsqcurvefit(@exponential, coeff0, time2, icelln_mcherry(i,:));
%             yhat_mcherry(i,:)=exponential(coeff_mcherry(i,:), time2);
%         else
%             continue
%         end
%     end
% end
% 
% figure,hold on
% for n=1:height(icelln_green)
%     if ismember(n, skip_green)==0
%         plot(time2, icelln_green(n,:))
%         %scatter(time2, yhat_green(n,:))
%     end
% end

%based on normalized 
% lin_green=[1 2 3 4 8 10 13 14];
% lin_mcherry=[];
% 
% exp_green=[5 6 7 9 11 12 15];
% exp_mcherry=[];
% 
% %define initial coeff
% coeff0a=[-1, 0.05];
% coeff0b=[1, 0.05];
% 
% %store fitresult
% expFit=cell(length(nGreen),1);
% expGOF=cell(length(nGreen),1);
% linFit=cell(length(nGreen),1);
% linGOF=cell(length(nGreen),1);
%     
% fitresult=cell(length(nGreen),1);
% gof=cell(length(nGreen),1);
%     
% for n=1:height(icelln_green)
%     
%     [xData, yData] = prepareCurveData(time2, icelln_green(n,:));
%     
%     % Set up fittype and options.
%     ft = fittype( 'exp1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.StartPoint = [1.01151960663985 -0.000417556351326276];
%     
%     % Fit model to data.
%     [fitresult{n,:}, gof{n,:}] = fit( xData, yData, ft, opts );
%     
%     % Plot fit with data.
%     figure( 'Name', num2str(n) );
%     h = plot( fitresult{n,:}, xData, yData );
%     legend( h, 'mNeonGreen Fluor Intesity (A.U.) vs. Time (s)', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
%     % Label axes
%     xlabel( 'time2', 'Interpreter', 'none' );
%     ylabel( 'icelln_green', 'Interpreter', 'none' );
%     grid on
%     pause, close
%     
% end

% for n=1:height(icelln_green)
%     if ismember(n,exp_green)
%         
%         [xData, yData] = prepareCurveData(time2, icelln_green(n,:));
%         
%         % Set up fittype and options.
%         ft = fittype( 'exp1' );
%         opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%         opts.Display = 'Off';
%         opts.StartPoint = [1.01151960663985 -0.000417556351326276];
%         
%         % Fit model to data.
%         [fitresult, gof] = fit( xData, yData, ft, opts );
%         
%         % Plot fit with data.
%         figure( 'Name', num2str(n) );
%         h = plot( fitresult, xData, yData );
%         legend( h, 'mNeonGreen Fluor Intesity (A.U.) vs. Time (s)', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
%         % Label axes
%         xlabel( 'time2', 'Interpreter', 'none' );
%         ylabel( 'icelln_green', 'Interpreter', 'none' );
%         grid on
%         pause, close
%         %[expFit{n,:}, linGOF{n,:}]=createExpFit(time2, icelln_green(n,:));
%         %yhat_green(n,:)=linear(coeff_green(n,:), time2);
%     elseif ismember(n, lin_green)
%         [xData, yData] = prepareCurveData(time2, icelln_green(n,:));
%         
%         % Set up fittype and options.
%         ft = fittype( 'poly1' );
% 
%         % Fit model to data.
%         [fitresult, gof] = fit( xData, yData, ft );
% 
%         % Plot fit with data.
%         figure( 'Name', num2str(n) );
%         h = plot( fitresult, xData, yData );
%         legend( h, 'mNeonGreen Fluor Intesity (A.U.) vs. Time (s)', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
%         
%         % Label axes
%         xlabel( 'time2', 'Interpreter', 'none' );
%         ylabel( 'icelln_green', 'Interpreter', 'none' );
%         grid on
%         pause, close
%         %[linFit{n,:}, linGOF{n,:}]=createLinFit(time2, icelln_green(n,:));
%         %yhat_green(n,:)=exponential(coeff_green(n,:), time2);
%     end
% end

% gof_green=[2,5,6,7,8,13,14];
% tconst=[];
% for n=1:height(fitresult)
%     if ismember(n, gof_green)==1
%         tconst=[tconst; fitresult{n,1}.b];
%     end
% end
% 
% figure
% histogram(tconst, 7)
% cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
% saveas(gcf, [basename '_histogram.fig'])
% saveas(gcf, [basename '_histogram.png'])

tconst=[];
for n=2:height(fitresult)
    tconst=[tconst; -1/fitresult{n,1}.b];
end

figure
histogram(tconst, height(icelln_green))
cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
xlabel('Time Scale Tau (s-)')
saveas(gcf, [basename '_histogram.fig'])
saveas(gcf, [basename '_histogram.png'])

figure, hold on
for n=2:height(icelln_green)
    plot(time2, icelln_green(n,:))
end
xlabel('time (s)')

ylabel('intensity (A.U.)')
cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
saveas(gcf, [basename '_histogramTraces.fig'])
saveas(gcf, [basename '_histogramTraces.png'])

%%%%%%%%%%%Functions
function [y] = exponential(b,x)
%this function calculates y=A*(e^-t/tau)
%where b(1)=A, b(2)=tau, x=t, and the cellular intensity=y;
y=b(1)*exp(-x.*b(2));
end

function [y] = linear(b,x)
%this function calculates y=A*(e^-t/tau)
%where b(1)=A, b(2)=tau, x=t, and the cellular intensity=y;
y=x.*b(1)+b(2);
end