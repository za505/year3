%Author: Zarina Akbary
%Date: 02/09/2022
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

% LBa=dir(['10232021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 1
% LBb=dir(['10262021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 2
% LB20=dir(['02042022_Exp1' '*dm.mat']); %LB, fr = 20 min
% LB5=dir(['02052022_Exp1' '*dm.mat']); %LB, fr = 5 min

LB1s=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate = 1.2 s
LB2s=dir(['12082021_Exp1' '*dm.mat']); %LB, frame rate = 2 s
LB3s=dir(['12082021_Exp2' '*dm.mat']); %LB, frame rate = 3 s

%% aggreate data
[adjintensity_LB1s, intensity_LB1s, lcell_LB1s, time_LB1s, tme_LB1s, imstart_LB1s]=controlAnalysis(LB1s);
[adjintensity_LB2s, intensity_LB2s, lcell_LB2s, time_LB2s, tme_LB2s, imstart_LB2s]=controlAnalysis(LB2s);
[adjintensity_LB3s, intensity_LB3s, lcell_LB3s, time_LB3s, tme_LB3s, imstart_LB3s]=controlAnalysis(LB3s);

%% fit background intensity to exponential
[bgintensity_LB1s, tau_LB1s, yhat_LB1s]=bgAnalysis(LB1s, imstart_LB1s);
[bgintensity_LB2s, tau_LB2s, yhat_LB2s]=bgAnalysis(LB2s, imstart_LB2s);
[bgintensity_LB3s, tau_LB3s, yhat_LB3s]=bgAnalysis(LB3s, imstart_LB3s);

%% plot 
figure, hold on
plot(tme_LB1s, (adjintensity_LB1s-bgintensity_LB1s(end)), '-r');
plot(tme_LB1s, bgintensity_LB1s, '-k');

% figure, hold on
% plot(tme_LB2s, adjintensity_LB2s, '-b');
% plot(tme_LB2s, bgintensity_LB2s, '-k');
% 
% figure, hold on
% plot(tme_LB3s, adjintensity_LB3s, '-g'); 
% plot(tme_LB3s, bgintensity_LB3s, '-k'); 

rstd_LB1s=stdFun(adjintensity_LB1s, bgintensity_LB1s);
rstd_LB2s=stdFun(adjintensity_LB2s, bgintensity_LB2s);
rstd_LB3s=stdFun(adjintensity_LB3s, bgintensity_LB3s);

% figure, hold on
% plot(tme_LB1s, yhat_LB1s, '--m');
% plot(tme_LB1s, bgintensity_LB1s, '-k');
% 
% figure, hold on
% plot(tme_LB2s, yhat_LB2s, '--m');
% plot(tme_LB2s, bgintensity_LB2s, '-k');
% 
% figure, hold on
% plot(tme_LB3s, yhat_LB3s, '--m'); 
% plot(tme_LB3s, bgintensity_LB3s, '-k'); 

dt1=tme_LB1s(end)-tme_LB1s(end-1);
tinterp=[tme_LB1s(end)+dt1:dt1:10];
interp_LB1s=interp1(tme_LB1s, bgintensity_LB1s, tinterp, 'linear', 'extrap');

% figure, hold on
% plot(tinterp,interp_LB1s, '--m');
% plot(tme_LB1s, bgintensity_LB1s, '-k');

x=tme_LB1s';
y=bgintensity_LB1s';
f = fit(x, y,'exp2');
yhat_LB1s=feval(f, 1:20);
figure, hold on
plot(1:20, yhat_LB1s, '--m');
plot(tme_LB1s, bgintensity_LB1s, '-k');

x=tme_LB2s';
y=bgintensity_LB2s';
f = fit(x, y,'exp2');
yhat_LB2s=feval(f, 1:20);
figure, hold on
plot(1:20, yhat_LB2s, '--m');
plot(tme_LB2s, bgintensity_LB2s, '-k');

x=tme_LB3s';
y=bgintensity_LB3s';
f = fit(x, y,'exp2');
yhat_LB3s=feval(f, 1:30);
figure, hold on
plot(1:100, yhat_LB3s, '--m');
plot(tme_LB3s, bgintensity_LB3s, '-k');

figure, hold on
plot(1:30, yhat_LB3s, '--m');
plot(tme_LB3s, adjintensity_LB3s, '-k');

x=tme_LB1s';
y=adjintensity_LB1s(1,:)';
f = fit(x, y,'exp1');
yhat_LB1s=feval(f, tme_LB1s);
figure, hold on
plot(tme_LB1s, yhat_LB1s, '--m');
plot(tme_LB1s, adjintensity_LB1s(1,:), '-k');
%% Functions
%rough standard deviation calculation
function rstd = stdFun(adjintensity, bgintensity)
    [nrow, ~]=size(adjintensity);
    A=(sum((adjintensity(:, end)-bgintensity(end)).^2))./nrow;
    rstd=sqrt(A);
end

%to aggregate control data and normalize to the post-lysis frame
function [adjintensity, intensity, lCell, time, tme, imstart]=controlAnalysis(datadir)
        
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
        
        %find the post-lysis frame 
        dt=round(diff(time), 2);
        if dt(1)==dt(end)
            dl=diff(lCell, 1, 2); 
            dlvg=mean(dl, 1, 'omitnan');
            [~, imstart] = min(dlvg);
            imstart=imstart+1; %this is at the detergent step, be careful 
        else
            if dt(end)<1 %discrepancies in dt might be found in the meta data, this is the ad hoc fix
                imstart=min(find(dt==dt(end)))+2;
            else 
                imstart=min(find(dt==dt(end)));
            end
        end
        
        %remove cells without a value for imstart
        idx=find(~isnan(intensity(:, imstart)));
        intensity=intensity(idx,:);
        
        %adjust the intensity values
        adjintensity=intensity(:, imstart:end);
            
        %adjust the time vector
        tme=tme(imstart:end)-tme(imstart);
                              
end

%to aggregate control background data and fit to an exponential
function [bgintensity, tau, yhat]=bgAnalysis(datadir, imstart)
        
        %pre-allocate variables
        bgintensity=[];
      
        %go through the data for each position
        for i=1:length(datadir)
            
            %load decayMeasure .mat file
            cd(datadir(i).folder)
            load(datadir(i).name, 'bg_intensity', 'time')

            bgintensity=[bgintensity; bg_intensity];

        end
        
        bgintensity=bgintensity(:, imstart:end);
        tme=time(imstart:end);

        modelfun = @(tau, x)(x(1)-x(end))*exp(-x./tau)+x(end);

        tau=nlinfit(tme, bgintensity, modelfun, 1);

        yhat=modelfun(tau, tme);
                              
end

%to aggregate data and normalize to the pre-lysis frame
function [normintensity, adjintensity, intensity, lCell, time, tme, imstart, imend]=preNormalize(datadir)
        
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
        intensity=intensity(idx,:);
        
        %subtract the final fluor. value
        adjintensity=intensity(:, imstart:end)-intensity(:, end);
                     
        %adjust the time vector
        tme=tme(imstart:end)-tme(imstart);
                         
        %normalize to the initial post-lysis frame
        normintensity=adjintensity./adjintensity(:,1);
        normintensity(normintensity<0)=NaN; 
        
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
        
        %find when the plateau starts
        imend=min(find(mean(adjintensity, 1, 'omitnan')<1));
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