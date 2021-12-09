%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all
%% matching fluorescent traces to images
%there are two questions I am trying to answer 1) do the cells with
%increasing fluorescent traces look different from cells with decreasing
%fluorescence? 2) what is the alpha that yields a straight line for the 1.2
%s frame rate LB control?

dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/12032021_analysis';
cd([dirsave '/MatFiles'])

datadir1=dir(['10232021_Exp1' '*dm.mat']); %LB, rep1
datadir2=dir(['10262021_Exp1' '*dm.mat']); %LB, rep2

datadir3=dir(['10232021_Exp2' '*dm.mat']); %PBS, 1 hour, rep 1
datadir4=dir(['10262021_Exp2' '*dm.mat']); %PBS, 1 hour, rep 2

datadir5=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate=1.2 s

%% calculate corrected fluorescence traces
%testing whether the correction is better when traces are normalized to
%initial post-lysis frame
[normintensity1, Cnew1, tme1, idx1, B1, positions1]=photoCorrect(datadir1, 10, 16.9441);
[normintensity2, Cnew2, tme2, idx2, B2, positoins2]=photoCorrect(datadir2, 9, 16.9441);

[normintensity3, Cnew3, tme3, idx3, B3, positions3]=photoCorrect(datadir3, 6, 16.9441);
[normintensity4, Cnew4, tme4, idx4, B4, positions4]=photoCorrect(datadir4, 9, 16.9441);

[normintensity5, Cnew5, tme5, idx5, B5, positions5]=photoCorrect(datadir5, 9, 16.9441);

%% plot normalized traces vs corrected traces
% comparePlot(normintensity1, Cnew1, tme1)
% comparePlot(normintensity2, Cnew2, tme2)
% 
% comparePlot(normintensity3, Cnew3, tme3)
% comparePlot(normintensity4, Cnew4, tme4)
% 
% comparePlot(normintensity5, Cnew5, tme5)
%% are the cells with increasing traces actually brighter? 
compareImage(datadir1, 10, idx1, B1, positions1, Cnew1)

%% Functions
function [normintensity, Cnew, tme, idx, bounds, positions]=photoCorrect(datadir, imstart, parameter)
        
        %pre-allocate variables
        intensity=[];
        bounds={};
        positions=nan(1, length(datadir));
        
        %go through the data for each position
        for i=1:length(datadir)
            
            %load decayMeasure .mat file
            cd(datadir(i).folder)
            load(datadir(i).name, 'icell_intensity', 'time', 'B')

            for n=1:height(icell_intensity)

                %make sure there are fluor readings during the initial frame and the final frame, otherwise the adjust. and norm. will look off
                if ~isnan(icell_intensity(n, imstart))&~isnan(icell_intensity(n, end)) 
                    intensity=[intensity; icell_intensity(n, :)];
                    bounds=[bounds; {B{n, imstart:end}}];
                end
            end

            if i==1
                tme=time(imstart:end)-time(imstart); %new time vector
            end
            
            positions(i) = height(intensity);
        end

        %adjust the background
        adjintensity = intensity-intensity(:, end);

        %normalize to the initial post-lysis frame
        normintensity=adjintensity(:, imstart:end)./adjintensity(:,imstart);
        normintensity(normintensity<0)=0;
        
        %fit normalized traces to exponential decay function
        modelfun=@(tau,x)exp(-x./tau);
        tau0=1;
        
        %pre-allocate variables
        fit=[];
        tau=[];
        yhat=[];

        for i=1:height(normintensity)
%             %this will skip the detergent perfusion reading
%             idx=find(normintensity(i,:)<=1); 
            tau_temp=nlinfit(tme, normintensity(i,:), modelfun, tau0);
            y_hat=modelfun(tau_temp, tme);   

%             fit=[fit, i];
            tau=[tau, tau_temp];
            yhat=[yhat; y_hat]; 

        end
        
        %calculate dt (there's variation of dt during initial values, so it's
        %easier to use end values)
        dt=tme(end)-tme(end-1);
        %beta = 1/((16.9441)*dt + 0.1823); %1/tau
        beta = 1/((parameter)*dt + 0.1823); %1/tau
        alpha=1/beta/dt; %tau*dt
        
        %assume that the initial 'measured' fluorescence values and corrected
        %fluor. values will be equal. I prefer to pre-allocate with nan in case 
        %some values are missing in the raw data
        Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
        Cnew(:, 1)=normintensity(:, 1);

        %pre-allocate variables
        dCB=nan(height(normintensity), length(tme)-1);
        dCT=nan(height(normintensity), length(tme)-1);
        dCP=nan(height(normintensity), length(tme)-1); %this is the dCP, or loss attributable to permeability

        unb_frac=nan(size(normintensity)); %fraction of unbleached fluor. 
        unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point

        Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
        Cbl_exp(:, 1)=0;

        mx=[];
        for n=1:height(normintensity)
           
            for i=1:length(tme)-1
                
                if normintensity(n,i)<0.05
                    continue
                end
                
                dCB(n,i) = normintensity(n,i)/(alpha); %this is the amount of photobleaching that occured in our measured value

                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value

                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
                
            end  
            
            %when does the loop break (Cnew no longer calculated)
            mx=[mx,i];
        end
        
        %mean time stamp for Cnew
        idx=round(mean(mx));
        
        gamma = nan(height(normintensity), 1);
    
        for n=1:height(normintensity)
            p = dCP(n,end)/dt;
            gamma(n,1) = p/Cnew(n, end-1);
        end
        
end

function comparePlot(normintensity, Cnew, tme)

    figure, hold on
    for n=1:height(normintensity)
        plot(tme, normintensity(n, :), '-b')
        plot(tme, Cnew(n, :), '-g')
    end
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')

end

function compareImage(datadir, imstart, idx, B, positions, Cnew)
      
    mdx=length(Cnew)/2;
    
    for i=1:length(datadir)
        
        if i==1
            s1=1;
            s2=positions(i);
        else
            s1=positions(i-1);
            s2=positions(i);
        end
        
        %first, load the decay measure data
        cd(datadir(i).folder)
        load(datadir(i).name, 'fluo_directory', 'channels')
        
        cd(channels{1});
        imagename=fluo_directory{1}(imstart + idx).name;
        im=imread(imagename);
        
        figure, imshow(im, []), hold on
        for n=s1:s2
            if Cnew(n, mdx)>1
                plot(B{n,mdx}(:,1),B{n,mdx}(:,2),'-r')
            else
                plot(B{n,mdx}(:,1),B{n,mdx}(:,2),'-b')
            end
        end
        
        pause, close all
    end
end