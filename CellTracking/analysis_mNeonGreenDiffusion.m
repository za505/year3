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

datadir1=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate = 1.2 s
datadir2=dir(['12082021_Exp1' '*dm.mat']); %LB, frame rate = 2 s
datadir3=dir(['12082021_Exp2' '*dm.mat']); %LB, frame rate = 3 s

%% calculate corrected fluorescence traces
%testing whether the correction is better when traces are normalized to
% %initial post-lysis frame
[normintensity1, adjintensity1, Cnew1, tme1, tau1, yhat1]=photoCorrect(datadir1, 9, 28.9210);
[normintensity2, adjintensity2, Cnew2, tme2, tau2, yhat2]=photoCorrect(datadir2, 9, 28.9210);
[normintensity3, adjintensity3, Cnew3, tme3, tau3, yhat3]=photoCorrect(datadir3, 9, 28.9210);

% %% plot normalized traces vs corrected traces
comparePlot(normintensity1, Cnew1, tme1)
comparePlot(normintensity2, Cnew2, tme2)
comparePlot(normintensity3, Cnew3, tme3)

comparePlot(normintensity1, yhat1, tme1)
comparePlot(normintensity2, yhat2, tme2)
comparePlot(normintensity3, yhat3, tme3)

% comparePlot(adjintensity1(1:9, 8:358), adjintensity2(:, 8:358), tme2)
%% calculate alpha based on tau values
% find the mean and standard deviation of the time constants
tau_means = [mean(tau1, 'omitnan'), mean(tau2, 'omitnan'), mean(tau3, 'omitnan')];
tau_std = [std(tau1, 0, 'omitnan'), std(tau2, 0, 'omitnan'), std(tau3, 0, 'omitnan')];

linearCoef1 = polyfit([repelem(1.2/60, length(tau1)), repelem(2/60, length(tau2)), repelem(3/60, length(tau3))],[tau1, tau2, tau3],1);
linearFit1= polyval(linearCoef1,[0, repelem(1.2/60, length(tau1)), repelem(2/60, length(tau2)), repelem(3/60, length(tau3))]);

cd(dirsave)
figure, hold on
scatter(repelem(1.2/60, length(tau1)), tau1, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(repelem(2/60, length(tau2)), tau2, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(repelem(3/60, length(tau3)), tau3, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')

scatter([1.2/60,2/60, 3/60], tau_means, 'MarkerFaceColor', 'black')
errorbar(1.2/60, tau_means(1), tau_std(1), 'Color', 'black')
errorbar(2/60, tau_means(2), tau_std(2), 'Color', 'black')
errorbar(3/60, tau_means(3), tau_std(3), 'Color', 'black')

plot([0, repelem(1.2/60, length(tau1)), repelem(2/60, length(tau2)), repelem(3/60, length(tau3))], linearFit1, '--b')
xlim([0, 0.06])
ylabel('\tau (min^{-1})')
xticks([0 1.2/60 2/60 3/60])
xticklabels({'0', '1.2 s', '2 s', '3 s'})
title('Tau vs Frame Rate')
% saveas(gcf, 'alpha.fig')
% saveas(gcf, 'alpha.png')


%% Functions
function [normintensity, adjintensity, Cnew, tme, tau, yhat]=photoCorrect(datadir, imstart, parameter)
        
        %pre-allocate variables
        intensity=[];
      
        %go through the data for each position
        for i=1:length(datadir)
            
            %load decayMeasure .mat file
            cd(datadir(i).folder)
            load(datadir(i).name, 'icell_intensity', 'time')

            for n=1:height(icell_intensity)

                %make sure there are fluor readings during the initial frame and the final frame, otherwise the adjust. and norm. will look off
                if ~isnan(icell_intensity(n, imstart))&~isnan(icell_intensity(n, end)) 
                    intensity=[intensity; icell_intensity(n, :)];
                end
            end

            if i==1
                tme=time(imstart:end)-time(imstart); %new time vector
            end

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