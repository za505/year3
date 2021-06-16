%decayMeasure.m
%Zarina Akbary, updated 05/25/21
%Calculates decay as a function of frame rate and exposure

clear, close all

%INSTRUCTIONS FOR USE:
%run photobleachMeasure.m first

%INPUT
%basename: experiments of interest
%dirname: where .mat files are stored

%OUTPUT:
%timescale=array of time scales
%Ts=array of T (number of frames)
%switchF=array of switchFrames
%A=array of A values (A*e^(-alpha*t)+y0)
%f=cell of coeff for exponential eqxn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basenames=["06062021_Exp1_mNeonGreen", "06062021_Exp1_mCherry"];
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/06062021_Exp1_colony1'];
recrunch=0;
vis=0;
B=length(basenames);%number of main directories to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(dirname)
    load(['06062021_dm.mat'])
else
    
%go to directory where .mat files are stored
cd(dirname)
load('06062021_Exp1_BTfluo.mat')

%pre-allocate variables
icell_green=[];
icell_mcherry=[];

%define new time variable to make things easier later
tpt=480;
cutoff=find(time==tpt);
tidx=tidx(tidx>=cutoff);
time=time(tidx)-tpt;

%split the cell intensity array into separate variables based on channel
for h=1:length(channels)
   for i=1:ncells
       if h==1 & sum(isnan(icell{1,h}(i,tidx)))<8 & mean(icell{1,h}(i,:),'omitnan')>10000
           icell_green=[icell_green; icell{1,h}(i,tidx)];
       elseif h==2 & sum(isnan(icell{1,h}(i,tidx)))<8 & mean(icell{1,h}(i,:),'omitnan')>10000
           icell_mcherry=[icell_mcherry; icell{1,h}(i,tidx)];
       end 
    end
end

%initialize coefficients
coeff0=[40000, 0.1];
    
%pre-allocate variables
coeff1=nan(height(icell_green), length(coeff0));
coeff2=nan(height(icell_mcherry), length(coeff0));
yhat1=nan(height(icell_green), length(time));
yhat2=nan(height(icell_mcherry), length(time));

%here, we can fit each cell intensity to a nonlinear exponential eqxn
for i=1:height(icell_green)
    coeff1(i,:)=nlinfit(time, icell_green(i,:), @exponential, coeff0);
    yhat1(i,:)=exponential(coeff1(i,:), time);
end
   
for i=1:height(icell_mcherry)
    coeff2(i,:)=nlinfit(time, icell_mcherry(i,:), @exponential, coeff0);
    yhat2(i,:)=exponential(coeff2(i,:), time);
end

%plot to see how well the eqxn fits the data
figure, hold on
for i=1:height(icell_green)
    plot(time, icell_green(i,:))
    scatter(time, yhat1(i,:))
end
saveas(gcf, [basename,'_expGreen.fig'])
saveas(gcf, [basename,'_expGreen.png'])
 
figure, hold on
for i=1:height(icell_mcherry)
    plot(time, icell_mcherry(i,:))
    scatter(time, yhat2(i,:))
end
saveas(gcf, [basename,'_expCherry.fig'])
saveas(gcf, [basename,'_expCherry.png'])

save(['dm'])

end

function [y] = exponential(b,x)
%this function calculates y=A*(e^-t/tau)
%where b(1)=A, b(2)=tau, x=t, and the cellular intensity=y;
y=b(1)*exp(-x./b(2));
end

