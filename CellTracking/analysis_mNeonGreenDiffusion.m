%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze decayMeasure (dm.mat) data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all
%% Generating Figure 1a
%to demonstrate that the cell wall is rate-limiting

% colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#E872D0'};
% colorcode2={'#A1D5F7', '#FFBB9E', '#FFE3A3', '#EAA6F7', '#CBF09C', '#C7EEFF', '#FFB3C1', '#FACAF0'};

colorcode={[0 0.45 0.74], [0.85 0.33 0.1], [0.49 0.18 0.56], [0.47 0.67 0.19], [0.47 0.67 0.19], [0.3 0.75 0.93], [0.64 0.08 0.18], [0.91 0.45 0.82]};
colorcode2={[0.63 0.84 0.97], [1 0.73 0.62], [0.92 0.65 0.97], [0.8 0.94 0.61], [0.8 0.94 0.61], [0.78 0.93 1], [1 0.7 0.76], [0.98 0.79 0.94]};

%first, direct the code to the location of the dm.mat files 
dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/12112021_analysis';
cd([dirsave '/MatFiles'])

LBa=dir(['10232021_Exp1' '*dm.mat']); %LB, rep 1
LBb=dir(['10262021_Exp1' '*dm.mat']); %LB, rep 2
PBS60a=dir(['10232021_Exp2' '*dm.mat']); %PBS 60 min, rep 1
PBS60b=dir(['10262021_Exp2' '*dm.mat']); %PBS 60 min, rep 2
PBS20a=dir(['11192021_Exp1' '*dm.mat']); %PBS 20 min, rep 1
PBS20b=dir(['11302021_Exp1' '*dm.mat']); %PBS 20 min, rep 2
PBS2a=dir(['11192021_Exp2' '*dm.mat']); %PBS 2 min, rep 1
PBS2b=dir(['12082021_Exp3' '*dm.mat']); %PBS 2 min, rep 2

% %% Perform analysis with traces normalized to post-lysis frame
%calculate corrected fluor traces
[normintensity_LBa, intensity_LBa, Cnew_LBa, tme_LBa, tau_LBa, yhat_LBa, dCB_LBa, dCT_LBa, dCP_LBa, Cbl_exp_LBa, unb_frac_LBa]=photoCorrect(LBa, 10, 28.9210);
[normintensity_LBb, intensity_LBb, Cnew_LBb, tme_LBb, tau_LBb, yhat_LBb, dCB_LBb, dCT_LBb, dCP_LBb, Cbl_exp_LBb, unb_frac_LBb]=photoCorrect(LBb, 8, 28.9210);
[normintensity_PBS60a, intensity_PBS60a, Cnew_PBS60a, tme_PBS60a, tau_PBS60a, yhat_PBS60a, dCB_PBS60a, dCT_PBS60a, dCP_PBS60a, Cbl_exp_PBS60a, unb_frac_PBS60a]=photoCorrect(PBS60a, 6, 28.9210);
[normintensity_PBS60b, intensity_PBS60b, Cnew_PBS60b, tme_PBS60b, tau_PBS60b, yhat_PBS60b, dCB_PBS60b, dCT_PBS60b, dCP_PBS60b, Cbl_exp_PBS60b, unb_frac_PBS60b]=photoCorrect(PBS60b, 7, 28.9210);
[normintensity_PBS20a, intensity_PBS20a, Cnew_PBS20a, tme_PBS20a, tau_PBS20a, yhat_PBS20a, dCB_PBS20a, dCT_PBS20a, dCP_PBS20a, Cbl_exp_PBS20a, unb_frac_PBS20a]=photoCorrect(PBS20a, 8, 28.9210);
[normintensity_PBS20b, intensity_PBS20b, Cnew_PBS20b, tme_PBS20b, tau_PBS20b, yhat_PBS20b, dCB_PBS20b, dCT_PBS20b, dCP_PBS20b, Cbl_exp_PBS20b, unb_frac_PBS20b]=photoCorrect(PBS20b, 7, 28.9210);
[normintensity_PBS2a, intensity_PBS2a, Cnew_PBS2a, tme_PBS2a, tau_PBS2a, yhat_PBS2a, dCB_PBS2a, dCT_PBS2a, dCP_PBS2a, Cbl_exp_PBS2a, unb_frac_PBS2a]=photoCorrect(PBS2a, 13, 28.9210);
[normintensity_PBS2b, intensity_PBS2b, Cnew_PBS2b, tme_PBS2b, tau_PBS2b, yhat_PBS2b, dCB_PBS2b, dCT_PBS2b, dCP_PBS2b, Cbl_exp_PBS2b, unb_frac_PBS2b]=photoCorrect(PBS2b, 13, 28.9210);
 
%% plot correction variables
figure, plot(tme_LBa(1:end-1), dCB_LBa)
figure, plot(tme_LBa(1:end-1), dCT_LBa)
figure, plot(tme_LBa(1:end-1), dCP_LBa)
figure, plot(tme_LBa, Cbl_exp_LBa)
figure, plot(tme_LBa, unb_frac_LBa)
% 
% figure, plot(tme_PBS20a(1:end-1), dCB_PBS20a)
% figure, plot(tme_PBS20a(1:end-1), dCT_PBS20a)
% figure, plot(tme_PBS20a(1:end-1), dCP_PBS20a)
% figure, plot(tme_PBS20a, Cbl_exp_PBS20a)
% figure, plot(tme_PBS20a, unb_frac_PBS20a)

% %plot normalized traces vs corrected traces
%cd([dirsave '/postLysis'])
% comparePlot(normintensity_LBa, Cnew_LBa, tme_LBa)
% %saveas(gcf, '10232021_Exp1_noAdjTraces.fig'), saveas(gcf, '10232021_Exp1_noAdjTraces.png'), close
% comparePlot(normintensity_LBb, Cnew_LBb, tme_LBb)
% %saveas(gcf, '10262021_Exp1_noAdjTraces.fig'), saveas(gcf, '10262021_Exp1_noAdjTraces.png'), close
% comparePlot(normintensity_PBS60a, Cnew_PBS60a, tme_PBS60a)
% %saveas(gcf, '10232021_Exp2_noAdjTraces.fig'), saveas(gcf, '10232021_Exp2_noAdjTraces.png'), close
% comparePlot(normintensity_PBS60b, Cnew_PBS60b, tme_PBS60b)
% %saveas(gcf, '10262021_Exp2_noAdjTraces.fig'), saveas(gcf, '10262021_Exp2_noAdjTraces.png'), close
% comparePlot(normintensity_PBS20a, Cnew_PBS20a, tme_PBS20a)
% %saveas(gcf, '11192021_Exp1_noAdjTraces.fig'), saveas(gcf, '11192021_Exp1_noAdjTraces.png'), close
% comparePlot(normintensity_PBS20b, Cnew_PBS20b, tme_PBS20b)
% %saveas(gcf, '11302021_Exp1_noAdjTraces.fig'), saveas(gcf, '11302021_Exp1_noAdjTraces.png'), close
% comparePlot(normintensity_PBS2a, Cnew_PBS2a, tme_PBS2a)
% %saveas(gcf, '11192021_Exp2_noAdjTraces.fig'), saveas(gcf, '11192021_Exp2_noAdjTraces.png'), close
% comparePlot(normintensity_PBS2b, Cnew_PBS2b, tme_PBS2b)
% %saveas(gcf, '12082021_Exp3_noAdjTraces.fig'), saveas(gcf, '12082021_Exp3_noAdjTraces.png'), close

%plot average raw fluor. intensities
cd(dirsave)

figure, hold on
ciplot(mean(normintensity_LBa, 1, 'omitnan')-std(normintensity_LBa, 0, 1, 'omitnan'), mean(normintensity_LBa, 1, 'omitnan')+std(normintensity_LBa, 0, 1, 'omitnan'), tme_LBa, colorcode2{1})
plot(tme_LBa, mean(normintensity_LBa, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(normintensity_LBb, 1, 'omitnan')-std(normintensity_LBb, 0, 1, 'omitnan'), mean(normintensity_LBb, 1, 'omitnan')+std(normintensity_LBb, 0, 1, 'omitnan'), tme_LBb, colorcode2{2})
plot(tme_LBb, mean(normintensity_LBb, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)

ciplot(mean(normintensity_PBS60a, 1, 'omitnan')-std(normintensity_PBS60a, 0, 1, 'omitnan'), mean(normintensity_PBS60a, 1, 'omitnan')+std(normintensity_PBS60a, 0, 1, 'omitnan'), tme_PBS60a, colorcode2{3})
plot(tme_PBS60a, mean(normintensity_PBS60a, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)

ciplot(mean(normintensity_PBS60b, 1, 'omitnan')-std(normintensity_PBS60b, 0, 1, 'omitnan'), mean(normintensity_PBS60b, 1, 'omitnan')+std(normintensity_PBS60b, 0, 1, 'omitnan'), tme_PBS60b, colorcode2{4})
plot(tme_PBS60b, mean(normintensity_PBS60b, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

ciplot(mean(normintensity_PBS20a, 1, 'omitnan')-std(normintensity_PBS20a, 0, 1, 'omitnan'), mean(normintensity_PBS20a, 1, 'omitnan')+std(normintensity_PBS20a, 0, 1, 'omitnan'), tme_PBS20a, colorcode2{5})
plot(tme_PBS20a, mean(normintensity_PBS20a, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(normintensity_PBS20b, 1, 'omitnan')-std(normintensity_PBS20b, 0, 1, 'omitnan'), mean(normintensity_PBS20b, 1, 'omitnan')+std(normintensity_PBS20b, 0, 1, 'omitnan'), tme_PBS20b, colorcode2{6})
plot(tme_PBS20b, mean(normintensity_PBS20b, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

ciplot(mean(normintensity_PBS2a, 1, 'omitnan')-std(normintensity_PBS2a, 0, 1, 'omitnan'), mean(normintensity_PBS2a, 1, 'omitnan')+std(normintensity_PBS2a, 0, 1, 'omitnan'), tme_PBS2a, colorcode2{7})
plot(tme_PBS2a, mean(normintensity_PBS2a, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)

ciplot(mean(normintensity_PBS2b, 1, 'omitnan')-std(normintensity_PBS2b, 0, 1, 'omitnan'), mean(normintensity_PBS2b, 1, 'omitnan')+std(normintensity_PBS2b, 0, 1, 'omitnan'), tme_PBS2b, colorcode2{8})
plot(tme_PBS2b, mean(normintensity_PBS2b, 1, 'omitnan'), 'Color', colorcode{8}, 'LineWidth', 1)

legend('LB, rep 1, std', 'LB, rep 1, mean', 'LB, rep 2, std', 'LB, rep 2, mean', 'PBS 1 hour, rep 1, std', 'PBS 1 hour, rep 1, mean', 'PBS 1 hour, rep 2, std', 'PBS 1 hour, rep 2, mean', 'PBS 20 min, rep 1, std', 'PBS 20 min, rep 1, mean', 'PBS 20 min, rep 2, std', 'PBS 20 min, rep 2, mean', 'PBS 2 min, rep 1, std', 'PBS 2 min, rep 1, mean', 'PBS 2 min, rep 2, std', 'PBS 2 min, rep 2, mean')
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
% saveas(gcf, 'rawIntensity.png')
% saveas(gcf, 'rawIntensity.fig')
%% Functions
function [time, lCell]=lengthView(datadir)
        %pre-allocate variables
        lCell=[];
      
        %go through the data for each position
        for i=1:length(datadir)
            
            %load decayMeasure .mat file
            cd(datadir(i).folder)
            load(datadir(i).name, 'lcell', 'time', 'icell_intensity', 'basename')

            for n=1:height(lcell)          
                lCell=[lCell; lcell(n, :)];
            end

        end
        
%         cd('../')
%         cd('./lengthTraces')
% 
%         figure, hold on
%         for i=1:height(lCell)
%             plot(time, lCell(i,:))
%         end
%         xlabel('Time (minutes)')
%         ylabel('Length \mum')
%         saveas(gcf, [basename '_lengthTraces.fig'])
%         saveas(gcf, [basename '_lengthTraces.png'])
%         pause(0.1), close

end

function [normintensity, intensity, Cnew, tme, tau, yhat, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(datadir, imstart, parameter)
        
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
        %adjintensity = intensity-intensity(:, end);

        %normalize to the initial post-lysis frame
        normintensity=intensity(:, imstart:end)./intensity(:,imstart);
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
                
%                 if normintensity(n,i)<0.05
%                     continue
%                 end
                
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