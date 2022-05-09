%Author: Zarina Akbary
%Date: 03/26/2022
%Purpose: To normalize diffusion data from dm.mat scripts

clear, close all

%directory guide
%untreated, 1 minute frame rate, 100% intensity = '10232021_Exp1' & '10262021_Exp1' 
%untreated, 1 minute frame rate, 20% intensity = '02212022_Exp2' 
%untreated, 5 minute frame rate, 100% intensity = '02122022_Exp1' 
%untreated, 10 minute frame rate, 100% intensity = '02122022_Exp2'
%untreated, 20 minute frame rate, 100% intensity = '02092022_Exp1'
%untreated, 20 minute frame rate, 20% intensity = '02212022_Exp1'
 
%2 minute PBS incubation, 1 minute frame rate, 100% intensity = '11192021_Exp2' & '12082021_Exp3' 
%20 minute PBS incubation, 1 minute frame rate, 100% intensity = '11192021_Exp1' & '11302021_Exp1'
%60 minute PBS incubation, 1 minute frame rate, 100% intensity = '10232021_Exp2' & '10262021_Exp2'
%120 minute PBS incubation, 1 minute frame rate, 100% intensity = '01142022_Exp1'

%LB + 20 mM Mg2+, 1 minute frame rate, 100% intensity = '01172022_Exp1' 
%LB + 10 mM EDTA, 1 minute frame rate, 100% intensity = '01172022_Exp2' 
%LB + 0.5 ug/mL tunicamycin, 1 minute frame rate, 100% intensity = '01242022_Exp1' 
%LB + 1 ug/mL vancomycin, 1 minute frame rate, 100% intensity = '01262022_Exp1'
%spent LB, 10 minute frame rate, 100% intensity = '02122022_Exp3' & '02192022_Exp1'
%spent LB, 1 minute frame rate, 100% intensity = '03012022_Exp1'

%untreated, 1.2 s frame rate, 100% intensity = '11202021_Exp1' 
%untreated, 2 s frame rate, 100% intensity = '12082021_Exp1'
%untreated, 3 s frame rate, 100% intensity = '12082021_Exp2'

%untreated, 1.76 s frame rate, 20% intensity = '02192022_Exp2' 
%untreated, 2.3 s frame rate, 20% intensity = '02192022_Exp3'
%untreated, 3 s frame rate, 20% intensity = '02192022_Exp4'

%prior to 03/26/2022
%for photobleach correction (100% intensity)
% alpha=32.2114;
% intercept=0.1614;

%for photobleach correction (20% intensity)
% alpha=151.7544;
% intercept=-1.0566;

%03/26/2022 analysis
%for photobleach correction (100% intensity)
% alpha=33.5358;
% intercept=-0.0204;

%for photobleach correction (20% intensity)
% alpha=140.7629;
% intercept=-0.8802;
%% User Input

%location and names of of the *.mat files 
dirsave='/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis';
basenames={'10232021_Exp1', '10262021_Exp1', '04042022_Exp1', '04052022_Exp3', '02212022_Exp2', '02122022_Exp1','02122022_Exp2', '02092022_Exp1', '02212022_Exp1', '11192021_Exp2', '12082021_Exp3', '11192021_Exp1', '11302021_Exp1', '10232021_Exp2', '10262021_Exp2', '01142022_Exp1', '01172022_Exp1', '01172022_Exp2', '01242022_Exp1', '01262022_Exp1', '02122022_Exp3', '02192022_Exp1', '03012022_Exp1', '04112022_Exp3', '11202021_Exp1', '12082021_Exp1', '12082021_Exp2', '04052022_Exp1', '04052022_Exp2', '02192022_Exp2', '02192022_Exp3', '02192022_Exp4'};

%% correct data
for i=1:length(basenames)
    
    basename=basenames{i};
    cd([dirsave '/normalizedFiles'])
    datadir=dir([basename '*']);
    
    power20 = {'02212022_Exp2', '02212022_Exp1', '02192022_Exp2', '02192022_Exp3', '02192022_Exp4'}; 
    
    if sum(strcmp(basename,power20))==1 % 20% intensity
    %for photobleach correction (20% intensity)
    alpha=140.7629;
    intercept=-0.8802;

%         alpha=151.7544;
%         intercept=-1.0566;
        
        load(datadir.name)
        [Cnew, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(tme, normintensity, alpha, intercept);

        cd([dirsave '/correctedFiles'])
        save([basename '_corrected.mat'], 'time', 'tme', 'alpha', 'intercept', 'Cnew', 'dCB', 'dCT', 'dCP', 'Cbl_exp', 'unb_frac');
    
    else
        
        %for photobleach correction (100% intensity)
        alpha=33.5358;
        %alpha=20;
        intercept=-0.0204;

%         alpha=32.2114;
%         intercept=0.1614;
        
        load(datadir.name);
        [Cnew, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(tme, normintensity, alpha, intercept);

        cd([dirsave '/correctedFiles'])
        save([basename '_corrected.mat'], 'time', 'tme', 'alpha', 'intercept', 'Cnew', 'dCB', 'dCT', 'dCP', 'Cbl_exp', 'unb_frac');
    end
    
end

%% Functions
%to correct for photobleaching
function [Cnew, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(tme, normintensity, alpha, intercept)
        
        %pre-allocate variables
        %assume that the initial 'measured' fluorescence values and corrected
        %fluor. values will be equal. I prefer to pre-allocate with nan in case 
        %some values are missing in the raw data
        Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
        Cnew(:, 1)=normintensity(:,1);
             
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
       %dC=@(C, alpha, dt, b)((C*exp(-dt/(alpha*dt+b)))*(1/(alpha*dt+b))*dt);
        
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