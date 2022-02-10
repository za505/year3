%Author: Zarina Akbary
%Date: 02/10/2022
%Purpose: To combine and analyze decayMeasure (dm.mat) data from mNeonGreen diffusion
%experiments. 

clear, close all

%Inputs
%dirsave = directory where all the dm.mat files are located and where plots
%are saved

%Outputs
%intensity = cell x time matrix of raw intensity values; concatenated
%icell_intensity variable from dm.mat files

%time = 1 x time matrix of time (minutes) values that correspond to
%intensity traces

%bgintensity = m x time matrix of background intensity values; concatenated
%bg_intensity variable from dm.mat files

%tme = 1 x n matrix of time (minutes) values that correspond to
%the pre- or post-lysis time point (imstart) and imend

%adjintensity = cell x tme maxtrix of intensity values minus bg intensity
%values from imstart to imend

%normintensity = cell x tme matrix of adjintensity values normalized to the
%initial value

%% User Input 

%color codes must be in RGB [0-1] format to be used in ciplot
colorcode={[204 0 0], [204 102 0], [204 204 0], [102 204 0], [0 204 204], [0 0 204], [102 0 204], [204 0 204], [204 0 102], [255 102 102], [255 178 102], [102 255 102], [102 255 255], [102 178 255], [178 102 255], [255 102 255],[255 102 178]};
colorcode2={[255 51 51], [255 153 51], [255 255 51], [153 255 51], [51 255 255], [51 51 255], [153 51 255], [255 51 255], [255 51 153], [255 204 204], [255 229 204], [204 255 204], [204 255 255], [204 229 255], [229 204 255], [255 204 255],[255 204 229]};

colorcode=cellfun(@(x)(x./255), colorcode, 'UniformOutput', false);
colorcode2=cellfun(@(x)(x./255), colorcode2, 'UniformOutput', false);

transparency = 0.3; %this is the alpha argument for the ciplot function

%location of the dm.mat files 
dirsave='/Users/zarina/Downloads/NYU/Year3_2022_Spring/02102022_reanalysis';
cd(dirsave)

LB1a=dir(['10232021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 1, full priming
LB1b=dir(['10262021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 2, full priming
LB1c=dir(['02092022_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 2, short priming
LB20a=dir(['02082022_Exp1' '*dm.mat']); %LB, fr = 20 min, full priming
LB20b=dir(['02042022_Exp1' '*dm.mat']); %LB, fr = 20 min, short priming

LB1s=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate = 1.2 s
LB2s=dir(['12082021_Exp1' '*dm.mat']); %LB, frame rate = 2 s
LB3s=dir(['12082021_Exp2' '*dm.mat']); %LB, frame rate = 3 s

%% calculate normalized fluorescence traces
[intensity_LB1a, bgintensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lcell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize1(LB1a, 1); %5
[intensity_LB1b, bgintensity_LB1b, adjintensity_LB1b, normintensity_LB1b, lcell_LB1b, time_LB1b, tme_LB1b, imstart_LB1b]=dataNormalize1(LB1b, 1); %4
% [intensity_LB1a, bgintensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lCell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize1(LB1c, 1); %5
% [intensity_LB1a, bgintensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lCell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize1(LB20a, 1); %5
% [intensity_LB1a, bgintensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lCell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize1(LB20b, 1); %5
% 
% [intensity_LB1a, bgintensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lCell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize1(LB1s, 2); 
% [intensity_LB1a, bgintensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lCell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize1(LB2s, 2); 
% [intensity_LB1a, bgintensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lCell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize1(LB3s, 2); 

%% combine datasets
normintensity_LB = [normintensity_LB1a(:, 1:length(tme_LB1b)); normintensity_LB1b];
adjintensity_LB = [adjintensity_LB1a(:, 1:length(time_LB1b)); adjintensity_LB1b];
bgintensity_LB = [bgintensity_LB1a(:, 1:length(time_LB1b)); bgintensity_LB1b];
intensity_LB = [intensity_LB1a(:, 1:length(time_LB1b)); intensity_LB1b];
lcell_LB = [lcell_LB1a(:, 1:length(time_LB1b)); lcell_LB1b];
time_LB=time_LB1b;
tme_LB=tme_LB1b;

%% plot length traces
figure, plot(time_LB, lcell_LB, '-b');
figure, plot(time_LB1c, lcell_LB1c, '-b');
figure, plot(time_LB20, lcell_LB20, '-b'); 

figure, hold on
ciplot(mean(lcell_LB, 1, 'omitnan')-std(lcell_LB, 0, 1, 'omitnan'), mean(lcell_LB, 1, 'omitnan')+std(lcell_LB, 0, 1, 'omitnan'), time_LB, colorcode2{1}, transparency)
plot(time_LB, mean(lcell_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(lcell_LB1c, 1, 'omitnan')-std(lcell_LB1c, 0, 1, 'omitnan'), mean(lcell_LB1c, 1, 'omitnan')+std(lcell_LB1c, 0, 1, 'omitnan'), time_LB1c, colorcode2{5}, transparency)
plot(time_LB1c, mean(lcell_LB1c, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(lcell_LB20, 1, 'omitnan')-std(lcell_LB20, 0, 1, 'omitnan'), mean(lcell_LB20, 1, 'omitnan')+std(lcell_LB20, 0, 1, 'omitnan'), time_LB20, colorcode2{3}, transparency)
plot(time_LB20, mean(lcell_LB20, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Length')

figure, hold on
plot(time_LB1s, lcell_LB1s, '-r');
plot(time_LB2s, lcell_LB2s, '-b');
plot(time_LB3s, lcell_LB3s, '-g'); 

%% plot the raw fluorescence intensity
figure, hold on, plot(time_LB, intensity_LB, '-b'), plot(time_LB, bgintensity_LB, '--k');
figure, plot(time_LB1c, intensity_LB1c, '-b');
figure, plot(time_LB20, intensity_LB20, '-b'); 

figure, hold on
ciplot(mean(intensity_LB, 1, 'omitnan')-std(intensity_LB, 0, 1, 'omitnan'), mean(intensity_LB, 1, 'omitnan')+std(intensity_LB, 0, 1, 'omitnan'), time_LB, colorcode2{1}, transparency)
plot(time_LB, mean(intensity_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(intensity_LB1c, 1, 'omitnan')-std(intensity_LB1c, 0, 1, 'omitnan'), mean(intensity_LB1c, 1, 'omitnan')+std(intensity_LB1c, 0, 1, 'omitnan'), time_LB1c, colorcode2{5}, transparency)
plot(time_LB1c, mean(intensity_LB1c, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(intensity_LB20, 1, 'omitnan')-std(intensity_LB20, 0, 1, 'omitnan'), mean(intensity_LB20, 1, 'omitnan')+std(intensity_LB20, 0, 1, 'omitnan'), time_LB20, colorcode2{3}, transparency)
plot(time_LB20, mean(intensity_LB20, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Length')

figure, hold on
plot(time_LB1s, intensity_LB1s, '-r');
plot(time_LB2s, intensity_LB2s, '-b');
plot(time_LB3s, intensity_LB3s, '-g'); 

%% plot adjusted intensity traces
figure, plot(time_LB, adjintensity_LB, '-r');
figure, plot(tme_LB1c, adjintensity_LB1c, '-r');
figure, plot(tme_LB20, adjintensity_LB20, '-r'); 

figure, hold on
plot(tme_LB1s, adjintensity_LB1s, '-r');
plot(tme_LB2s, adjintensity_LB2s, '-b');
plot(tme_LB3s, adjintensity_LB3s, '-g'); 

%% plot normalized fluorescence traces
figure, plot(tme_LB, normintensity_LB, '-g');
figure, plot(tme_LB1c, normintensity_LB1c, '-g');
figure, plot(tme_LB20, normintensity_LB20, '-g'); 

figure, hold on
ciplot(mean(normintensity_LB, 1, 'omitnan')-std(normintensity_LB, 0, 1, 'omitnan'), mean(normintensity_LB, 1, 'omitnan')+std(normintensity_LB, 0, 1, 'omitnan'), tme_LB, colorcode2{1}, transparency)
plot(tme_LB, mean(normintensity_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(normintensity_LB1c, 1, 'omitnan')-std(normintensity_LB1c, 0, 1, 'omitnan'), mean(normintensity_LB1c, 1, 'omitnan')+std(normintensity_LB1c, 0, 1, 'omitnan'), tme_LB1c, colorcode2{5}, transparency)
plot(tme_LB1c, mean(normintensity_LB1c, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(normintensity_LB20, 1, 'omitnan')-std(normintensity_LB20, 0, 1, 'omitnan'), mean(normintensity_LB20, 1, 'omitnan')+std(normintensity_LB20, 0, 1, 'omitnan'), tme_LB20, colorcode2{3}, transparency)
plot(tme_LB20, mean(normintensity_LB20, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)

xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_LB1s, normintensity_LB1s, '-r');
plot(tme_LB2s, normintensity_LB2s, '-b');
plot(tme_LB3s, normintensity_LB3s, '-g'); 

%% calculate alpha
[tau1, yhat_LB1s]=tauCalc(tme_LB1s, normintensity_LB1s);
[tau2, yhat_LB2s]=tauCalc(tme_LB2s, normintensity_LB2s);
[tau3, yhat_LB3s]=tauCalc(tme_LB3s, normintensity_LB3s);

figure, hold on
plot(tme_LB1s, normintensity_LB1s, '-k')
plot(tme_LB1s, yhat_LB1s, '--r')

figure, hold on
plot(tme_LB2s, normintensity_LB2s, '-k')
plot(tme_LB2s, yhat_LB2s, '--r')

figure, hold on
plot(tme_LB3s, normintensity_LB3s, '-k')
plot(tme_LB3s, yhat_LB3s, '--r')

tau_means = [mean(tau1, 'omitnan'), mean(tau2, 'omitnan'), mean(tau3, 'omitnan')];
tau_std = [std(tau1, 0, 'omitnan'), std(tau2, 0, 'omitnan'), std(tau3, 0, 'omitnan')];

dt1=tme_LB1s(end)-tme_LB1s(end-1);
dt2=tme_LB2s(end)-tme_LB2s(end-1);
dt3=tme_LB3s(end)-tme_LB3s(end-1);

linearCoef1 = polyfit([repelem(dt1, length(tau1)), repelem(dt2, length(tau2)), repelem(dt3, length(tau3))],[tau1', tau2', tau3'],1);
linearFit1= polyval(linearCoef1,[0 dt1 dt2 dt3]);

figure, hold on
scatter(repelem(dt1, length(tau1)), tau1, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(repelem(dt2, length(tau2)), tau2, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(repelem(dt3, length(tau3)), tau3, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')

scatter([dt1, dt2, dt3], tau_means, 'MarkerFaceColor', 'black')
errorbar(dt1, tau_means(1), tau_std(1), 'Color', 'black')
errorbar(dt2, tau_means(2), tau_std(2), 'Color', 'black')
errorbar(dt3, tau_means(3), tau_std(3), 'Color', 'black')

plot([0 dt1 dt2 dt3], linearFit1, '--b')
xlim([0, 0.06])
ylabel('\tau (min^{-1})')
xticks([0 dt1 dt2 dt3])
xticklabels({'0', '1.2 s', '2 s', '3 s'})
title('Tau vs Frame Rate')
 
%% correct for photobleaching
% alpha=38.898;
% intercept=-0.0765;

% alpha=28.9210;
% intercept=0.2482;

alpha=39.1226;
intercept=-0.0827;

[Cnew_LB, dCB_LB, dCT_LB, dCP_LB, CblExp_LB, unbFrac_LB]=photoCorrect(tme_LB, normintensity_LB, alpha, intercept);
[Cnew_LB1c, dCB_LB1c, dCT_LB1c, dCP_LB1c, CblExp_LB1c, unbFrac_LB1c]=photoCorrect(tme_LB1c, normintensity_LB1c, alpha, intercept);
[Cnew_LB20, dCB_LB20, dCT_LB20, dCP_LB20, CblExp_LB20, unbFrac_LB20]=photoCorrect(tme_LB20, normintensity_LB20, alpha, intercept);

[Cnew_LB1s, dCB_LB1s, dCT_LB1s, dCP_LB1s, CblExp_LB1s, unbFrac_LB1s]=photoCorrect(tme_LB1s, normintensity_LB1s, alpha, intercept);
[Cnew_LB2s, dCB_LB2s, dCT_LB2s, dCP_LB2s, CblExp_LB2s, unbFrac_LB2s]=photoCorrect(tme_LB2s, normintensity_LB2s, alpha, intercept);
[Cnew_LB3s, dCB_LB3s, dCT_LB3s, dCP_LB3s, CblExp_LB3s, unbFrac_LB3s]=photoCorrect(tme_LB3s, normintensity_LB3s, alpha, intercept);

%% plot corrected traces
figure, plot(tme_LB, Cnew_LB, '-r');
figure, plot(tme_LB1c, Cnew_LB1c, '-r');
figure, plot(tme_LB20, Cnew_LB20, '-r'); 

figure, hold on
ciplot(mean(Cnew_LB, 1, 'omitnan')-std(Cnew_LB, 0, 1, 'omitnan'), mean(Cnew_LB, 1, 'omitnan')+std(Cnew_LB, 0, 1, 'omitnan'), tme_LB, colorcode2{1}, transparency)
plot(tme_LB, mean(Cnew_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(Cnew_LB1c, 1, 'omitnan')-std(Cnew_LB1c, 0, 1, 'omitnan'), mean(Cnew_LB1c, 1, 'omitnan')+std(Cnew_LB1c, 0, 1, 'omitnan'), tme_LB1c, colorcode2{5}, transparency)
plot(tme_LB1c, mean(Cnew_LB1c, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(Cnew_LB20, 1, 'omitnan')-std(Cnew_LB20, 0, 1, 'omitnan'), mean(Cnew_LB20, 1, 'omitnan')+std(Cnew_LB20, 0, 1, 'omitnan'), tme_LB20, colorcode2{3}, transparency)
plot(tme_LB20, mean(Cnew_LB20, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_LB1s, Cnew_LB1s, '-r');
plot(tme_LB2s, Cnew_LB2s, '-b');
plot(tme_LB3s, Cnew_LB3s, '-g'); 

%% Functions
%to aggregate data, subtract background, and normalize
function [intensity, bgintensity, adjintensity, normintensity, lCell, time, tme, imstart]=dataNormalize1(datadir, set)
        
        %set 1=pre-lysis
        %set 2=post-lysis 
        
        %pre-allocate variables
        intensity=[];
        lCell=[];
        bgintensity=[];
        adjintensity=[];
        
        %go through the data for each position
        for i=1:length(datadir)
            
            %load decayMeasure .mat file
            cd(datadir(i).folder)
            load(datadir(i).name, 'icell_intensity', 'bg_intensity', 'time', 'lcell')
            
            adj_intensity=icell_intensity-bg_intensity;
            adjintensity=[adjintensity; adj_intensity];
            intensity=[intensity; icell_intensity];
            lCell=[lCell; lcell];
            bgintensity=[bgintensity; bg_intensity];
            
            if i==1
                tme=time; %pre-set new time vector
            end

        end
        
        if set==1
            %find the final pre-lysis frame 
            dl=diff(lCell, 1, 2);
            lvg=mean(dl, 1, 'omitnan');
            [~, imstart]=min(lvg);
        else
            %find the initial post-lysis frame 
            dt=round(diff(time), 2);
            if dt(1)==dt(end)
                dl=diff(lCell, 1, 2); 
                dlvg=mean(dl, 1, 'omitnan');
                [~, imstart] = min(dlvg);
                imstart=imstart+2;
            else
                if dt(end)<1 %discrepancies in dt might be found in the meta data, this is the ad hoc fix
                    imstart=min(find(dt==dt(end)))+2;
                else 
                    imstart=min(find(dt==dt(end)));
                end
            end
        end
        
        %remove cells without a value for imstart
        idx=find(~isnan(intensity(:, imstart)));
        intensity=intensity(idx,:);
        adjintensity=adjintensity(idx, :);
        lCell=lCell(idx,:);
            
        %set the limit of detection at 200 AU
        %adjintensity(adjintensity<200)=NaN;
        %adjintensity(adjintensity<0)=NaN;

        %find when all values are below the limit of detection
%         [nrow, ncol]=size(adjintensity);
%         nsum=sum(isnan(adjintensity), 1);
%         if max(find(nsum==0))>1
%             imend=max(find(nsum==0));
%         else 
%             imend=ncol;
%         end
        
%         %adjust the time points in the fluor matrix
%         adjintensity=adjintensity(:, imstart:imend);

        %adjust the time vector
        tme=tme(imstart:end)-tme(imstart);
                         
        %normalize to the initial post-lysis frame
        normintensity=adjintensity(:, imstart:end)./adjintensity(:,1);              
        
end

%to aggregate data and normalize to the pre-lysis frame
function [normintensity, adjintensity, intensity, lCell, time, tme, imstart, imend]=dataNormalize2(datadir)
        
        %pre-allocate variables
        intensity=[];
        lCell=[];
      
        %go through the data for each position
        for i=1:length(datadir)
            
            %load decayMeasure .mat file
            cd(datadir(i).folder)
            load(datadir(i).name, 'icell_intensity', 'time', 'lcell')

            intensity=[intensity; icell_intensity];
            lCell=[lCell; lcell];
     
            if i==1
                tme=time; %pre-set new time vector
            end

        end
        
        %find the final pre-lysis frame 
        dl=diff(lCell, 1, 2);
        lvg=mean(dl, 1, 'omitnan');
        [~, imstart]=min(lvg);
        
        %remove cells without a value for imstart
        idx=find(~isnan(intensity(:, imstart)));
        adjintensity=intensity(idx,:);
        
        %subtract the final fluor. value
        adjintensity=adjintensity-adjintensity(:, end);
        adjintensity(adjintensity<0)=NaN; 

        %find when none of the values are below the limit of detection
        [nrow, ncol]=size(adjintensity);
        nsum=sum(isnan(adjintensity), 1);
        if max(find(nsum==0))>1
            imend=max(find(nsum==0));
        else 
            imend=ncol;
        end
        
        %adjust the time points in the fluor matrix
        adjintensity=adjintensity(:, imstart:imend);

        %adjust the time vector
        tme=tme(imstart:imend)-tme(imstart);
                         
        %normalize to the initial pre-lysis frame
        normintensity=adjintensity./adjintensity(:,1);  
        
        %interpolate the fluor values during detergent perfusion
        for n=1:height(normintensity)
            idx=find(normintensity(n,:)>1);
            if ~isempty(idx)
                [nrow, ncol]=size(normintensity);
                x=setdiff(1:ncol, idx); %x=time
                v=normintensity(n, x); %v=intensity
                vq=interp1(x, v, idx); %vq=interpolation at query time pts
                normintensity(n,idx)=vq;
            end
        end
        
%         %find when the plateau starts
%         imend=min(find(mean(adjintensity, 1, 'omitnan')<1));
end

%to correct for photobleaching
function [Cnew, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(tme, normintensity, alpha, intercept)
        
        %pre-allocate variables
        %assume that the initial 'measured' fluorescence values and corrected
        %fluor. values will be equal. I prefer to pre-allocate with nan in case 
        %some values are missing in the raw data
        Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
        Cnew(:, 1)=normintensity(:, 1);
        
        
        dCB=nan(height(normintensity), length(tme)-1); %change in fluor. due to photobleaching
        dCT=nan(height(normintensity), length(tme)-1); %total change in fluor.
        dCP=nan(height(normintensity), length(tme)-1); %this is the dCP, or loss attributable to permeability

        unb_frac=nan(size(normintensity)); %fraction of unbleached fluor. 
        unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point

        Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
        Cbl_exp(:, 1)=0;
        
        %calculate dt (the dt between frames may vary in a single run).
        %Note: this is also the frame rate
        dt=round(diff(tme), 2);
        
        %this formula comes from the slope and intercept calculated for the 1.2, 2,
        %and 3 second tau vs frame rate controls 
       dC=@(C, alpha, dt, b)(C/(alpha*dt+b))*dt;
       %dC=@(C, alpha, intercept, dt, b)(C/(alpha*dt+b));
        
        %the correction
        for n=1:height(normintensity)
           
            for i=1:length(tme)-1
                
                dCB(n,i) = dC(normintensity(n,i), alpha, dt(i), intercept); %this is the amount of photobleaching that occured in our measured value

                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value

                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)*unb_frac(n,i); %Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
                
            end  
            
        end      

end

%do not fit to y=(1-beta)*exp(-t./tau)+beta (beta=normalized limit of detection)
%fit to y=Ae^(-t/tau)
function [tau, yhat]=tauCalc(tme, normintensity)
    
    [nrow, ~]=size(normintensity);
    tau=nan(nrow, 1);
    modelfun = @(tau, x)exp(-x./tau);
    for i=1:nrow
        tau(i, 1)=nlinfit(tme, normintensity(i,:), modelfun, 1);
    end
    
    yhat=modelfun(tau, tme);
    
end