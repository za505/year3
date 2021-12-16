%Author: Zarina Akbary
%Date: 12/15/2021
%Purpose: To generate and analyze data for Figure 1
%to demonstrate that the cell wall is rate-limiting

clear, close all

%INSTRUCTIONS FOR USE:
%run BacTrack.m first

%INPUT
%basename: experiment of interest
%dirname: where BT.mat files are stored
%dirsave: where the .mat output of this script is saved

%OUTPUT:
%B=BacTrack.m variable, boundary of tracked cells
%T=BacTrack.m variable, total time points
%pixels=BacTrack.m variable, coordinates of pixels that correspond to
%tracked cells
%time=BacTrack.m variable, time vector
%lcell=BacTrack.m variable, length traces of tracked cells
%imstart=initial post-lysis frame
%tme=indexed 'time' variable from imstart to end
%cell_intensity=fluorescence intensity of tracked cells
%bg_intensity=pixel values of selected background
%% User Input
basename=['12082021_Exp3'];
fluordir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12082021_analysis/12082021_Exp3/' basename '_colony3/' basename '_mNeonGreen/' basename '_aligned'];
fluorsave=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12082021_analysis/12082021_Exp3/' basename '_colony3/' basename '_mNeonGreen/' basename '_figures'];
phasedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12082021_analysis/12082021_Exp3/' basename '_colony3/' basename '_phase/' basename '_figures'];
dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/12152021_analysis';
imstart=13;
recrunch=0;

% colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#E872D0'};
% colorcode2={'#A1D5F7', '#FFBB9E', '#FFE3A3', '#EAA6F7', '#CBF09C', '#C7EEFF', '#FFB3C1', '#FACAF0'};

colorcode={[0 0.45 0.74], [0.85 0.33 0.1], [0.49 0.18 0.56], [0.47 0.67 0.19], [0.47 0.67 0.19], [0.3 0.75 0.93], [0.64 0.08 0.18], [0.91 0.45 0.82]};
colorcode2={[0.63 0.84 0.97], [1 0.73 0.62], [0.92 0.65 0.97], [0.8 0.94 0.61], [0.8 0.94 0.61], [0.78 0.93 1], [1 0.7 0.76], [0.98 0.79 0.94]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    load([basename '_colony3_dm.mat'])
else

    %load fluorescent images
    cd(fluordir); 
    fluo_directory=dir('*.tif');

    %load BacTrack data
    cd(phasedir)
    load([basename '_BTphase.mat'], 'B', 'T', 'time', 'pixels', 'ncells')
    time=time/60;
    
    %preallocate cells
    cellIntensity=nan(ncells,T); %cellular intensities
    %cellAdjIntensity=nan(ncells,T); %cellular intensities-bg
    %cellNormIntensity=nan(ncells, T-imstart);
    bgIntensity=nan(1,T); %background intensity

    %determine region where you'll measure background intensity
    cd(fluordir)
    imagename=fluo_directory(1).name;
    im=imread(imagename);    
    [p1, p2]=getBackground(imagename);
    close all

    %now let's track intensity over time
    for t=1:T

        t

        %load the image
        imagename=fluo_directory(t).name;
        im=imread(imagename);

        %measure background level
        bglevel = measureBackground(imagename, p1, p2);
        bgIntensity(t)=bglevel; 

        %measure cellular intensity
        for j=1:height(pixels)

            %calculate intensity
            cellIntensity(j,t)=mean(im(pixels{j,t}));

        end

    end

    %substract the background intensity
    cellAdjIntensity=cellIntensity-bgIntensity;

    %normalize fluor. based on initial post-lysis intensity
    cellNormIntensity=cellAdjIntensity(:, imstart:end)./cellAdjIntensity(:, imstart);
    tme=time(imstart:end)-time(imstart);
    
    %calculate average intensity and standard deviation
    avgIntensity = mean(cellIntensity, 1, 'omitnan');
    stdIntensity = std(cellIntensity, 0, 1, 'omitnan');

    %calculate average adjusted intensity and standard deviation
    avgAdjIntensity = mean(cellAdjIntensity, 1, 'omitnan');
    stdAdjIntensity = std(cellAdjIntensity, 0, 1, 'omitnan');

    %calculate average normalized intensity and standard deviation
    avgNormIntensity = mean(cellNormIntensity, 1, 'omitnan');
    stdNormIntensity = std(cellNormIntensity, 0, 1, 'omitnan');
end
   
%save the variables
cd(dirsave)
save([basename '_colony3_dm.mat'])

%Plot data
%cd(fluorsave)

figure, hold on
for n=1:ncells
    plot(time, cellIntensity(n, :), 'Color', colorcode{1})
end
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
ylim([0 Inf])
saveas(gcf, [basename '_colony3_cellIntensity.png'])
saveas(gcf, [basename '_colony3_cellIntensity.fig'])

figure, hold on
for n=1:ncells
    plot(time, cellAdjIntensity(n, :), 'Color', colorcode{2})
end
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence (A.U.)')
ylim([0 Inf])
saveas(gcf, [basename '_colony3_cellAdjIntensity.png'])
saveas(gcf, [basename '_colony3_cellAdjIntensity.fig'])

figure, hold on
for n=1:ncells
    plot(tme, cellNormIntensity(n, :), 'Color', colorcode{5})
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
ylim([0 Inf])
saveas(gcf, [basename '_colony3_cellNormIntensity.png'])
saveas(gcf, [basename '_colony3_cellNormIntensity.fig'])

figure
plot(time, bgIntensity, 'Color', colorcode{3})
xlabel('Time (minutes)')
ylabel('Background Fluorescence (A.U.)')
ylim([0 Inf])
saveas(gcf, [basename '_colony3_bgIntensity.png'])
saveas(gcf, [basename '_colony3_bgIntensity.fig'])

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

function bglevel = measureBackground(imagename, p1, p2)
        
        %Load last image
        %imagename=fluo_directory{i}(t).name;
        im2=imread(imagename);
        
        %Determine background
        backim=im2(p1(2):p2(2),p1(1):p2(1));
%         [counts,bins]=imhist(backim);
%         [~,binnum]=max(counts);
%         maxpos=bins(binnum);
        bglevel=mean(mean(backim));
        
end 

%% Code Graveyard
% %first, direct the code to the location of the dm.mat files 
% dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/12112021_analysis';
% cd([dirsave '/MatFiles'])
% 
% LBa=dir(['12082021_Exp3' '*dm.mat']); %LB, rep 1
% LBb=dir(['12082021_Exp3' '*dm.mat']); %LB, rep 2
% PBS60a=dir(['12082021_Exp3' '*dm.mat']); %PBS 60 min, rep 1
% PBS60b=dir(['12082021_Exp3' '*dm.mat']); %PBS 60 min, rep 2
% PBS20a=dir(['12082021_Exp3' '*dm.mat']); %PBS 20 min, rep 1
% PBS20b=dir(['12082021_Exp3' '*dm.mat']); %PBS 20 min, rep 2
% PBS2a=dir(['12082021_Exp3' '*dm.mat']); %PBS 2 min, rep 1
% PBS2b=dir(['12082021_Exp3' '*dm.mat']); %PBS 2 min, rep 2
% 
% %plot average raw fluor. intensities
% % cd(dirsave)
% % 
% % figure, hold on
% % ciplot(mean(intensity_LBa, 1, 'omitnan')-std(intensity_LBa, 0, 1, 'omitnan'), mean(intensity_LBa, 1, 'omitnan')+std(intensity_LBa, 0, 1, 'omitnan'), time_LBa, colorcode2{1})
% % plot(time_LBa, mean(intensity_LBa, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
% % 
% % ciplot(mean(intensity_LBb, 1, 'omitnan')-std(intensity_LBb, 0, 1, 'omitnan'), mean(intensity_LBb, 1, 'omitnan')+std(intensity_LBb, 0, 1, 'omitnan'), time_LBb, colorcode2{2})
% % plot(time_LBb, mean(intensity_LBb, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)
% % 
% % ciplot(mean(intensity_PBS60a, 1, 'omitnan')-std(intensity_PBS60a, 0, 1, 'omitnan'), mean(intensity_PBS60a, 1, 'omitnan')+std(intensity_PBS60a, 0, 1, 'omitnan'), time_PBS60a, colorcode2{3})
% % plot(time_PBS60a, mean(intensity_PBS60a, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
% % 
% % ciplot(mean(intensity_PBS60b, 1, 'omitnan')-std(intensity_PBS60b, 0, 1, 'omitnan'), mean(intensity_PBS60b, 1, 'omitnan')+std(intensity_PBS60b, 0, 1, 'omitnan'), time_PBS60b, colorcode2{4})
% % plot(time_PBS60b, mean(intensity_PBS60b, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)
% % 
% % ciplot(mean(intensity_PBS20a, 1, 'omitnan')-std(intensity_PBS20a, 0, 1, 'omitnan'), mean(intensity_PBS20a, 1, 'omitnan')+std(intensity_PBS20a, 0, 1, 'omitnan'), time_PBS20a, colorcode2{5})
% % plot(time_PBS20a, mean(intensity_PBS20a, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)
% % 
% % ciplot(mean(intensity_PBS20b, 1, 'omitnan')-std(intensity_PBS20b, 0, 1, 'omitnan'), mean(intensity_PBS20b, 1, 'omitnan')+std(intensity_PBS20b, 0, 1, 'omitnan'), time_PBS20b, colorcode2{6})
% % plot(time_PBS20b, mean(intensity_PBS20b, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
% % 
% % ciplot(mean(intensity_PBS2a, 1, 'omitnan')-std(intensity_PBS2a, 0, 1, 'omitnan'), mean(intensity_PBS2a, 1, 'omitnan')+std(intensity_PBS2a, 0, 1, 'omitnan'), time_PBS2a, colorcode2{7})
% % plot(time_PBS2a, mean(intensity_PBS2a, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)
% % 
% % ciplot(mean(intensity_PBS2b, 1, 'omitnan')-std(intensity_PBS2b, 0, 1, 'omitnan'), mean(intensity_PBS2b, 1, 'omitnan')+std(intensity_PBS2b, 0, 1, 'omitnan'), time_PBS2b, colorcode2{8})
% % plot(time_PBS2b, mean(intensity_PBS2b, 1, 'omitnan'), 'Color', colorcode{8}, 'LineWidth', 1)
% % 
% % legend('LB, rep 1, std', 'LB, rep 1, mean', 'LB, rep 2, std', 'LB, rep 2, mean', 'PBS 1 hour, rep 1, std', 'PBS 1 hour, rep 1, mean', 'PBS 1 hour, rep 2, std', 'PBS 1 hour, rep 2, mean', 'PBS 20 min, rep 1, std', 'PBS 20 min, rep 1, mean', 'PBS 20 min, rep 2, std', 'PBS 20 min, rep 2, mean', 'PBS 2 min, rep 1, std', 'PBS 2 min, rep 1, mean', 'PBS 2 min, rep 2, std', 'PBS 2 min, rep 2, mean')
% % xlabel('Time (minutes)')
% % ylabel('Fluorescence (A.U.)')
% % saveas(gcf, 'rawIntensity.png')
% % saveas(gcf, 'rawIntensity.fig')
% 
% 
% % %% Perform analysis with traces normalized to post-lysis frame
% %calculate corrected fluor traces
% [normintensity_LBa, intensity_LBa, Cnew_LBa, tme_LBa, tau_LBa, yhat_LBa, dCB_LBa, dCT_LBa, dCP_LBa, Cbl_exp_LBa, unb_frac_LBa]=photoCorrect(LBa, 10, 28.9210);
% [normintensity_LBb, intensity_LBb, Cnew_LBb, tme_LBb, tau_LBb, yhat_LBb, dCB_LBb, dCT_LBb, dCP_LBb, Cbl_exp_LBb, unb_frac_LBb]=photoCorrect(LBb, 8, 28.9210);
% [normintensity_PBS60a, intensity_PBS60a, Cnew_PBS60a, tme_PBS60a, tau_PBS60a, yhat_PBS60a, dCB_PBS60a, dCT_PBS60a, dCP_PBS60a, Cbl_exp_PBS60a, unb_frac_PBS60a]=photoCorrect(PBS60a, 6, 28.9210);
% [normintensity_PBS60b, intensity_PBS60b, Cnew_PBS60b, tme_PBS60b, tau_PBS60b, yhat_PBS60b, dCB_PBS60b, dCT_PBS60b, dCP_PBS60b, Cbl_exp_PBS60b, unb_frac_PBS60b]=photoCorrect(PBS60b, 7, 28.9210);
% [normintensity_PBS20a, intensity_PBS20a, Cnew_PBS20a, tme_PBS20a, tau_PBS20a, yhat_PBS20a, dCB_PBS20a, dCT_PBS20a, dCP_PBS20a, Cbl_exp_PBS20a, unb_frac_PBS20a]=photoCorrect(PBS20a, 8, 28.9210);
% [normintensity_PBS20b, intensity_PBS20b, Cnew_PBS20b, tme_PBS20b, tau_PBS20b, yhat_PBS20b, dCB_PBS20b, dCT_PBS20b, dCP_PBS20b, Cbl_exp_PBS20b, unb_frac_PBS20b]=photoCorrect(PBS20b, 7, 28.9210);
% [normintensity_PBS2a, intensity_PBS2a, Cnew_PBS2a, tme_PBS2a, tau_PBS2a, yhat_PBS2a, dCB_PBS2a, dCT_PBS2a, dCP_PBS2a, Cbl_exp_PBS2a, unb_frac_PBS2a]=photoCorrect(PBS2a, 13, 28.9210);
% [normintensity_PBS2b, intensity_PBS2b, Cnew_PBS2b, tme_PBS2b, tau_PBS2b, yhat_PBS2b, dCB_PBS2b, dCT_PBS2b, dCP_PBS2b, Cbl_exp_PBS2b, unb_frac_PBS2b]=photoCorrect(PBS2b, 13, 28.9210);
%  
% %% plot correction variables
% figure, plot(tme_LBa(1:end-1), dCB_LBa)
% figure, plot(tme_LBa(1:end-1), dCT_LBa)
% figure, plot(tme_LBa(1:end-1), dCP_LBa)
% figure, plot(tme_LBa, Cbl_exp_LBa)
% figure, plot(tme_LBa, unb_frac_LBa)
% % 
% % figure, plot(tme_PBS20a(1:end-1), dCB_PBS20a)
% % figure, plot(tme_PBS20a(1:end-1), dCT_PBS20a)
% % figure, plot(tme_PBS20a(1:end-1), dCP_PBS20a)
% % figure, plot(tme_PBS20a, Cbl_exp_PBS20a)
% % figure, plot(tme_PBS20a, unb_frac_PBS20a)
% 
% % %plot normalized traces vs corrected traces
% % cd([dirsave '/postLysis'])
% % comparePlot(intensity_LBa, Cnew_LBa, tme_LBa)
% % saveas(gcf, '12082021_Exp3_rawCorrectedTraces.fig'), saveas(gcf, '12082021_Exp3_rawCorrectedTraces.png'), close
% % comparePlot(normintensity_LBb, Cnew_LBb, tme_LBb)
% % saveas(gcf, '12082021_Exp3_rawCorrectedTraces.fig'), saveas(gcf, '12082021_Exp3_rawCorrectedTraces.png'), close
% % comparePlot(normintensity_PBS60a, Cnew_PBS60a, tme_PBS60a)
% % saveas(gcf, '12082021_Exp3_rawCorrectedTraces.fig'), saveas(gcf, '12082021_Exp3_rawCorrectedTraces.png'), close
% % comparePlot(normintensity_PBS60b, Cnew_PBS60b, tme_PBS60b)
% % saveas(gcf, '12082021_Exp3_rawCorrectedTraces.fig'), saveas(gcf, '12082021_Exp3_rawCorrectedTraces.png'), close
% % comparePlot(normintensity_PBS20a, Cnew_PBS20a, tme_PBS20a)
% % saveas(gcf, '12082021_Exp3_rawCorrectedTraces.fig'), saveas(gcf, '12082021_Exp3_rawCorrectedTraces.png'), close
% % comparePlot(normintensity_PBS20b, Cnew_PBS20b, tme_PBS20b)
% % saveas(gcf, '12082021_Exp3_rawCorrectedTraces.fig'), saveas(gcf, '12082021_Exp3_rawCorrectedTraces.png'), close
% % comparePlot(normintensity_PBS2a, Cnew_PBS2a, tme_PBS2a)
% % saveas(gcf, '12082021_Exp3_rawCorrectedTraces.fig'), saveas(gcf, '12082021_Exp3_rawCorrectedTraces.png'), close
% % comparePlot(normintensity_PBS2b, Cnew_PBS2b, tme_PBS2b)
% % saveas(gcf, '12082021_Exp3_rawCorrectedTraces.fig'), saveas(gcf, '12082021_Exp3_rawCorrectedTraces.png'), close
% 
% %% Functions
% function [time, lCell]=lengthView(datadir)
%         %pre-allocate variables
%         lCell=[];
%       
%         %go through the data for each position
%         for i=1:length(datadir)
%             
%             %load decayMeasure .mat file
%             cd(datadir(i).folder)
%             load(datadir(i).name, 'lcell', 'time', 'icell_intensity', 'basename')
% 
%             for n=1:height(lcell)          
%                 lCell=[lCell; lcell(n, :)];
%             end
% 
%         end
%         
% %         cd('../')
% %         cd('./lengthTraces')
% % 
% %         figure, hold on
% %         for i=1:height(lCell)
% %             plot(time, lCell(i,:))
% %         end
% %         xlabel('Time (minutes)')
% %         ylabel('Length \mum')
% %         saveas(gcf, [basename '_lengthTraces.fig'])
% %         saveas(gcf, [basename '_lengthTraces.png'])
% %         pause(0.1), close
% 
% end
% 
% function [normintensity, intensity, Cnew, tme, tau, yhat, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(datadir, imstart, parameter)
%         
%         %pre-allocate variables
%         intensity=[];
%         tau=[];
%         yhat=[];
%         
%         %go through the data for each position
%         for i=1:length(datadir)
%             
%             %load decayMeasure .mat file
%             cd(datadir(i).folder)
%             load(datadir(i).name, 'icell_intensity', 'time')
% 
%             for n=1:height(icell_intensity)
% 
%                 %make sure there are fluor readings during the initial frame and the final frame, otherwise the adjust. and norm. will look off
%                 if ~isnan(icell_intensity(n, imstart))&~isnan(icell_intensity(n, end)) 
%                     intensity=[intensity; icell_intensity(n, :)];
%                 end
%             end
% 
% %             if i==1
% %                 tme=time(imstart:end)-time(imstart); %new time vector
% %             end
%             
%             if i==1
%                 tme=time; %new time vector
%             end
% 
%         end
% 
%         normintensity=intensity;
% %         %adjust the background
% %         adjintensity = intensity-intensity(:, end);
% % 
% %         %normalize to the initial post-lysis frame
% %         normintensity=adjintensity(:, imstart:end)./adjintensity(:,imstart);
% %         normintensity(normintensity<0)=0;
% %         
% %         %fit normalized traces to exponential decay function
% %         modelfun=@(tau,x)exp(-x./tau);
% %         tau0=1;
% %         
% %         %pre-allocate variables
% %         fit=[];
% %         tau=[];
% %         yhat=[];
% % 
% %         for i=1:height(normintensity)
% % %             %this will skip the detergent perfusion reading
% % %             idx=find(normintensity(i,:)<=1); 
% %             tau_temp=nlinfit(tme, normintensity(i,:), modelfun, tau0);
% %             y_hat=modelfun(tau_temp, tme);   
% % 
% % %             fit=[fit, i];
% %             tau=[tau, tau_temp];
% %             yhat=[yhat; y_hat]; 
% % 
% %         end
% %     
%         %calculate dt (there's variation of dt during initial values, so it's
%         %easier to use end values)
%         dt=tme(end)-tme(end-1);
%         %beta = 1/((16.9441)*dt + 0.1823); %1/tau
%         beta = 1/((parameter)*dt + 0.2482); %1/tau
%         alpha=1/beta/dt; %tau*dt
%         
%         %assume that the initial 'measured' fluorescence values and corrected
%         %fluor. values will be equal. I prefer to pre-allocate with nan in case 
%         %some values are missing in the raw data
%         Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
%         Cnew(:, 1)=normintensity(:, 1);
% 
%         %pre-allocate variables
%         dCB=nan(height(normintensity), length(tme)-1);
%         dCT=nan(height(normintensity), length(tme)-1);
%         dCP=nan(height(normintensity), length(tme)-1); %this is the dCP, or loss attributable to permeability
% 
%         unb_frac=nan(size(normintensity)); %fraction of unbleached fluor. 
%         unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point
% 
%         Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
%         Cbl_exp(:, 1)=0;
% 
%         for n=1:height(normintensity)
%            
%             for i=1:length(tme)-1
%                 
% %                 if normintensity(n,i)<0.05
% %                     continue
% %                 end
%                 
%                 dCB(n,i) = normintensity(n,i)/(alpha); %this is the amount of photobleaching that occured in our measured value
% 
%                 dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
% 
%                 dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability
% 
%                 dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached
% 
%                 Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);
% 
%                 Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
%                 
%                 unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
%                 
%             end  
%             
%         end
%         
% %         gamma = nan(height(normintensity), 1);
% %     
% %         for n=1:height(normintensity)
% %             p = dCP(n,end)/dt;
% %             gamma(n,1) = p/Cnew(n, end-1);
% %         end
%         
% end
% 
% function comparePlot(normintensity, Cnew, tme)
% 
%     figure, hold on
%     for n=1:height(normintensity)
%         plot(tme, normintensity(n, :), '-b')
%         plot(tme, Cnew(n, :), '-g')
%     end
%     xlabel('Time (minutes)')
%     ylabel('Fluorescence (A.U.)')
% 
% end
% 
% function compareImage(datadir, imstart, idx, B, positions, Cnew)
%       
%     mdx=length(Cnew)/2;
%     
%     for i=1:length(datadir)
%         
%         if i==1
%             s1=1;
%             s2=positions(i);
%         else
%             s1=positions(i-1);
%             s2=positions(i);
%         end
%         
%         %first, load the decay measure data
%         cd(datadir(i).folder)
%         load(datadir(i).name, 'fluo_directory', 'channels')
%         
%         cd(channels{1});
%         imagename=fluo_directory{1}(imstart + idx).name;
%         im=imread(imagename);
%         
%         figure, imshow(im, []), hold on
%         for n=s1:s2
%             if Cnew(n, mdx)>1
%                 plot(B{n,mdx}(:,1),B{n,mdx}(:,2),'-r')
%             else
%                 plot(B{n,mdx}(:,1),B{n,mdx}(:,2),'-b')
%             end
%         end
%         
%         pause, close all
%     end
% end