%Author: Zarina Akbary
%Date: 02/07/2022
%Purpose: To combine and analyze decayMeasure (dm.mat) data from mNeonGreen diffusion
%experiments. 

clear, close all

%Inputs
%dirsave = directory where all the dm.mat files are located and where plots
%are saved

%Outputs
%intensity = m x n matrix of raw intensity values; concatenated
%icell_intensity variable from dm.mat files
%time = 1 x n matrix of time (minutes) values that correspond to
%intensity traces

%normintensity = m x n matrix of intensity values where 1) the background
%was subtracted 2) the intensity values during detergent perfusion were
%interpolated and 3) the intensity values were normalized by the initial
%pre-lysis value
%time = 1 x n matrix of time (minutes) values that correspond to
%normalized intensity traces
%Cnew = m x n matrix of normalized intensity values corrected for
%photobleaching
%beta = m x 1 matrix of final fluorescent values for corrected traces
%dCB = m x n matrix of change in corrected fluorescence attributable to
%photobleaching between two adjacent time points
%dCT = m x n matrix of change in corrected fluorescence between two time points
%dCP = m x n matrix of change in corrected fluorescence attributable to
%permeability between two adjacent time points
%CblExp = m x n matrix of the concentration of bleached fluorophores over
%time
%unbFrac = m x n matrix of the fraction of unbleached fluor over time
%midx = index of Cnew traces that correspond to normintensity traces; indexed
%Cnew traces have a beta value <=0.9 (greater than 0.9 are outliers with
%noisy corrected traces) and are not NaN
%tau = m x 1 matrix of time constants calculated by fitting Cnew to the
%exponential decay function y = (1-beta)*e^-t/tau+beta
%yhat = m x n matrix of predicted trace value
%initial = m x n matrix of initial phase intensities
%final = m x n matrix of final phase intensities
%gamma = m x 1 matrix of final/initial phase intensity ratios

%% User Input 

%color codes must be in RGB [0-1] format to be used in ciplot
colorcode={[204 0 0], [204 102 0], [204 204 0], [102 204 0], [0 204 204], [0 0 204], [102 0 204], [204 0 204], [204 0 102], [255 102 102], [255 178 102], [102 255 102], [102 255 255], [102 178 255], [178 102 255], [255 102 255],[255 102 178]};
colorcode2={[255 51 51], [255 153 51], [255 255 51], [153 255 51], [51 255 255], [51 51 255], [153 51 255], [255 51 255], [255 51 153], [255 204 204], [255 229 204], [204 255 204], [204 255 255], [204 229 255], [229 204 255], [255 204 255],[255 204 229]};

colorcode=cellfun(@(x)(x./255), colorcode, 'UniformOutput', false);
colorcode2=cellfun(@(x)(x./255), colorcode2, 'UniformOutput', false);

transparency = 0.3; %this is the alpha argument for the ciplot function

%location of the dm.mat files 
dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/02072022_analysis';
cd([dirsave '/MatFiles'])

LBa=dir(['10232021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 1
LBb=dir(['10262021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 2
LB20=dir(['02042022_Exp1' '*dm.mat']); %LB, fr = 20 min
LB5=dir(['02052022_Exp1' '*dm.mat']); %LB, fr = 5 min

LB1s=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate = 1.2 s
LB2s=dir(['12082021_Exp1' '*dm.mat']); %LB, frame rate = 2 s
LB3s=dir(['12082021_Exp2' '*dm.mat']); %LB, frame rate = 3 s

%% calculate normalized fluorescence traces
[normintensity_LBa, intensity_LBa, lcell_LBa, time_LBa, tme_LBa, imstart_LBa, imend_LBa]=dataNormalize(LBa); %5
[normintensity_LBb, intensity_LBb, lcell_LBb, time_LBb, tme_LBb, imstart_LBb, imend_LBb]=dataNormalize(LBb); %4
[normintensity_LB20, intensity_LB20, lcell_LB20, time_LB20, tme_LB20, imstart_LB20, imend_LB20]=dataNormalize(LB20); %5
[normintensity_LB5, intensity_LB5, lcell_LB5, time_LB5, tme_LB5, imstart_LB5, imend_LB5]=dataNormalize(LB5); %5

[normintensity_LB1s, intensity_LB1s, lcell_LB1s, time_LB1s, tme_LB1s, imstart_LB1s, imend_LB1s]=dataNormalize(LB1s); 
[normintensity_LB2s, intensity_LB2s, lcell_LB2s, time_LB2s, tme_LB2s, imstart_LB2s, imend_LB2s]=dataNormalize(LB2s); 
[normintensity_LB3s, intensity_LB3s, lcell_LB3s, time_LB3s, tme_LB3s, imstart_LB3s, imend_LB3s]=dataNormalize(LB3s); 

%% combine datasets
normintensity_LB = [normintensity_LBa(:, 1:length(tme_LBb)); normintensity_LBb];
intensity_LB = [intensity_LBa(:, 1:length(time_LBb)); intensity_LBb];
lcell_LB = [lcell_LBa(:, 1:length(time_LBb)); lcell_LBb];
time_LB=time_LBb;
tme_LB=tme_LBb;

%% plot length traces
figure, plot(time_LB, lcell_LB, '-b');
figure, plot(time_LB5, lcell_LB5, '-b');
figure, plot(time_LB20, lcell_LB20, '-b'); 

figure, hold on
ciplot(mean(lcell_LB, 1, 'omitnan')-std(lcell_LB, 0, 1, 'omitnan'), mean(lcell_LB, 1, 'omitnan')+std(lcell_LB, 0, 1, 'omitnan'), time_LB, colorcode2{1}, transparency)
plot(time_LB, mean(lcell_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(lcell_LB5, 1, 'omitnan')-std(lcell_LB5, 0, 1, 'omitnan'), mean(lcell_LB5, 1, 'omitnan')+std(lcell_LB5, 0, 1, 'omitnan'), time_LB5, colorcode2{5}, transparency)
plot(time_LB5, mean(lcell_LB5, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(lcell_LB20, 1, 'omitnan')-std(lcell_LB20, 0, 1, 'omitnan'), mean(lcell_LB20, 1, 'omitnan')+std(lcell_LB20, 0, 1, 'omitnan'), time_LB20, colorcode2{3}, transparency)
plot(time_LB20, mean(lcell_LB20, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Length')

figure, hold on
plot(time_LB1s, lcell_LB1s, '-r');
plot(time_LB2s, lcell_LB2s, '-b');
plot(time_LB3s, lcell_LB3s, '-g'); 

%% plot the raw fluorescence intensity
figure, plot(time_LB, intensity_LB, '-b');
figure, plot(time_LB5, intensity_LB5, '-b');
figure, plot(time_LB20, intensity_LB20, '-b'); 

figure, hold on
ciplot(mean(intensity_LB, 1, 'omitnan')-std(intensity_LB, 0, 1, 'omitnan'), mean(intensity_LB, 1, 'omitnan')+std(intensity_LB, 0, 1, 'omitnan'), time_LB, colorcode2{1}, transparency)
plot(time_LB, mean(intensity_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(intensity_LB5, 1, 'omitnan')-std(intensity_LB5, 0, 1, 'omitnan'), mean(intensity_LB5, 1, 'omitnan')+std(intensity_LB5, 0, 1, 'omitnan'), time_LB5, colorcode2{5}, transparency)
plot(time_LB5, mean(intensity_LB5, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(intensity_LB20, 1, 'omitnan')-std(intensity_LB20, 0, 1, 'omitnan'), mean(intensity_LB20, 1, 'omitnan')+std(intensity_LB20, 0, 1, 'omitnan'), time_LB20, colorcode2{3}, transparency)
plot(time_LB20, mean(intensity_LB20, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Length')

figure, hold on
plot(time_LB1s, intensity_LB1s, '-r');
plot(time_LB2s, intensity_LB2s, '-b');
plot(time_LB3s, intensity_LB3s, '-g'); 

%% plot normalized fluorescence traces
figure, plot(tme_LB, normintensity_LB, '-g');
figure, plot(tme_LB5, normintensity_LB5, '-g');
figure, plot(tme_LB20, normintensity_LB20, '-g'); 


figure, hold on
ciplot(mean(normintensity_LB, 1, 'omitnan')-std(normintensity_LB, 0, 1, 'omitnan'), mean(normintensity_LB, 1, 'omitnan')+std(normintensity_LB, 0, 1, 'omitnan'), tme_LB, colorcode2{1}, transparency)
plot(tme_LB, mean(normintensity_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(normintensity_LB5, 1, 'omitnan')-std(normintensity_LB5, 0, 1, 'omitnan'), mean(normintensity_LB5, 1, 'omitnan')+std(normintensity_LB5, 0, 1, 'omitnan'), tme_LB5, colorcode2{5}, transparency)
plot(tme_LB5, mean(normintensity_LB5, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(normintensity_LB20, 1, 'omitnan')-std(normintensity_LB20, 0, 1, 'omitnan'), mean(normintensity_LB20, 1, 'omitnan')+std(normintensity_LB20, 0, 1, 'omitnan'), tme_LB20, colorcode2{3}, transparency)
plot(tme_LB20, mean(normintensity_LB20, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)

xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_LB1s, normintensity_LB1s, '-r');
plot(tme_LB2s, normintensity_LB2s, '-b');
plot(tme_LB3s, normintensity_LB3s, '-g'); 

%% correct for photobleaching
alpha=28.9210;
[Cnew_LB, dCB_LB, dCT_LB, dCP_LB, CblExp_LB, unbFrac_LB]=photoCorrect(tme_LB, normintensity_LB, alpha);
[Cnew_LB5, dCB_LB5, dCT_LB5, dCP_LB5, CblExp_LB5, unbFrac_LB5]=photoCorrect(tme_LB5, normintensity_LB5, alpha);
[Cnew_LB20, dCB_LB20, dCT_LB20, dCP_LB20, CblExp_LB20, unbFrac_LB20]=photoCorrect(tme_LB20, normintensity_LB20, alpha);

[Cnew_LB1s, dCB_LB1s, dCT_LB1s, dCP_LB1s, CblExp_LB1s, unbFrac_LB1s]=photoCorrect(tme_LB1s, normintensity_LB1s, alpha);
[Cnew_LB2s, dCB_LB2s, dCT_LB2s, dCP_LB2s, CblExp_LB2s, unbFrac_LB2s]=photoCorrect(tme_LB2s, normintensity_LB2s, alpha);
[Cnew_LB3s, dCB_LB3s, dCT_LB3s, dCP_LB3s, CblExp_LB3s, unbFrac_LB3s]=photoCorrect(tme_LB3s, normintensity_LB3s, alpha);

%% plot corrected traces
figure, plot(tme_LB, Cnew_LB, '-r');
figure, plot(tme_LB5, Cnew_LB5, '-r');
figure, plot(tme_LB20, Cnew_LB20, '-r'); 


figure, hold on
ciplot(mean(Cnew_LB, 1, 'omitnan')-std(Cnew_LB, 0, 1, 'omitnan'), mean(Cnew_LB, 1, 'omitnan')+std(Cnew_LB, 0, 1, 'omitnan'), tme_LB, colorcode2{1}, transparency)
plot(tme_LB, mean(Cnew_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(Cnew_LB5, 1, 'omitnan')-std(Cnew_LB5, 0, 1, 'omitnan'), mean(Cnew_LB5, 1, 'omitnan')+std(Cnew_LB5, 0, 1, 'omitnan'), tme_LB5, colorcode2{5}, transparency)
plot(tme_LB5, mean(Cnew_LB5, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(Cnew_LB20, 1, 'omitnan')-std(Cnew_LB20, 0, 1, 'omitnan'), mean(Cnew_LB20, 1, 'omitnan')+std(Cnew_LB20, 0, 1, 'omitnan'), tme_LB20, colorcode2{3}, transparency)
plot(tme_LB20, mean(Cnew_LB20, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_LB1s, Cnew_LB1s, '-r');
plot(tme_LB2s, Cnew_LB2s, '-b');
plot(tme_LB3s, Cnew_LB3s, '-g'); 

%% Functions
%to normalize experimental traces with consistent frame rates
function [normintensity, intensity, lCell, time, tme, imstart, imend]=dataNormalize(datadir)
        
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
        
        %find the pre-lysis frame (the max length)
        dl=diff(lCell, 1, 2); 
        dlvg=mean(dl, 1, 'omitnan');
        [~, imstart] = min(dlvg);
        imstart=imstart-1;
        
        %remove cells without a value for imstart
        idx=find(~isnan(intensity(:, imstart)));
        intensity=intensity(idx,:);
        
        %set the limit of detection ignore values below it
        lmt=1400;
        adjintensity=intensity;
        adjintensity(adjintensity<=lmt)=-1; %to distinguish between true NaN values and those below the limit of detection 
        
        %find the end time point (where most of the cells are below the
        %limit of detection)
        nsum=sum(adjintensity, 1);
        if min(nsum) < 0
            [~, imend]=min(nsum);
        else
            [~, imend]=size(adjintensity);
        end
            
        %adjust the time vector
        tme=tme(imstart:imend)-tme(imstart);
                         
        %normalize to the initial pre-lysis frame
        normintensity=adjintensity(:, imstart:imend)./adjintensity(:,imstart);
        normintensity(normintensity<0)=NaN; %there are some cells that may have reached the limit of detection sooner than others
        
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
        
end

%to normalize control traces without consistent frame rates
function [normintensity, intensity, lCell, time, tme, imstart, imend]=dataNormalize(datadir)
        
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
        
        %find the pre-lysis frame (the max length)
        dl=diff(lCell, 1, 2); 
        dlvg=mean(dl, 1, 'omitnan');
        [~, imstart] = min(dlvg);
        imstart=imstart-1;
        
        %remove cells without a value for imstart
        idx=find(~isnan(intensity(:, imstart)));
        intensity=intensity(idx,:);
        
        %set the limit of detection ignore values below it
        lmt=1400;
        adjintensity=intensity;
        adjintensity(adjintensity<=lmt)=-1; %to distinguish between true NaN values and those below the limit of detection 
        
        %find the end time point (where most of the cells are below the
        %limit of detection)
        nsum=sum(adjintensity, 1);
        if min(nsum) < 0
            [~, imend]=min(nsum);
        else
            [~, imend]=size(adjintensity);
        end
            
        %adjust the time vector
        tme=tme(imstart:imend)-tme(imstart);
                         
        %normalize to the initial pre-lysis frame
        normintensity=adjintensity(:, imstart:imend)./adjintensity(:,imstart);
        normintensity(normintensity<0)=NaN; %there are some cells that may have reached the limit of detection sooner than others
        
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
        
end

%to correct for photobleaching
function [Cnew, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(tme, normintensity, alpha)
        
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
        dt=diff(tme);
        
        %this formula comes from the slope and intercept calculated for the 1.2, 2,
        %and 3 second tau vs frame rate controls 
       dC=@(C, alpha, dt, b)(C/(alpha*dt+b))*dt;
        
        %the correction
        for n=1:height(normintensity)
           
            for i=1:length(tme)-1
                
                dCB(n,i) = dC(normintensity(n,i), alpha, dt(i), 0.2482); %this is the amount of photobleaching that occured in our measured value

                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value

                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)*unb_frac(n,i); %Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
                
            end  
            
        end      

end

