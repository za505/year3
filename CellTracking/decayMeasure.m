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
basenames=["05272021_FITCK_e100", "05272021_FITCK_e001", "05272021_FITCK_e010"];
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05272021_analysis/05272021_FITCK_figures'];
recrunch=0;
vis=0;
B=length(basenames);%number of main directories to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%go to directory where .mat files are stored
cd(dirname)

%pre-allocate variables
timescales=nan(B,1);
Ts=nan(B,1);
time=cell(B,1);
intensity=cell(B,1);
coeff=nan(B,3);
yhat=cell(B,1);
A=nan(B,1);
f=cell(B,1);
yhat2=cell(B,1);
thalf=nan(B,1);
mindiff=nan(B,1);
exp2=[];

%initialize coeff
coeff0=[40000, -0.1, 0.1];

%let's go through the mat files and store the variables we need
for b=1:B
    base=char(basenames(b))
    load([base '_pb'])

    timescales(b,1)=tscale;
    Ts(b,1)=T;
    times{b,1}=time;
    intensity{b,1}=intensityAvg;
    A(b,1)=intensityAvg(1,1);
    
    %here, we can fit to a nonlinear exponential eqxn
    coeff(b,:)=nlinfit(times{b,1}, intensity{b,1}, @exponential, coeff0);
    
    %calculate the predicted y values
    yhat{b,1}=exponential(coeff(b,:), times{b,1});
    
    %plot to see how well the eqxn fits the data
    figure
    scatter(times{b,1}, intensity{b,1}), hold on
    plot(times{b,1}, yhat{b,1})
    pause
    saveas(gcf, [basename,'_exp1.fig'])
    saveas(gcf, [basename,'_exp1.png'])
    close
    
    %plot to see how well the eqxn fits the data
    f{b,1}=fit(times{b,1}', intensity{b,1}', 'exp2');
    figure
    plot(f{b,1}, times{b,1}, intensity{b,1})
    pause
    saveas(gcf, [basename,'_exp2.fig'])
    saveas(gcf, [basename,'_exp2.png'])
    close
    
    %calculate new yhat values
    yhat2{b,1}=f{b,1}(times{b,1});
    
    %pull out the coefficient values
    exp2(b,:)=coeffvalues(f{b,1});
end

%plot exposure as a function of intensity
exposure=[0.1, 1, 10];
k1=exp2(:,2);
k2=exp2(:,4);

yyaxis left
plot(exposure, k1)
ylabel('k1'), hold on
yyaxis right
plot(exposure,k2)
ylabel('k2')
title('Coefficients as a function of Exposure')

save([basename '_dm'])

function [y] = exponential(b,x)
%this function calculates y=A*(e^alpha*t)+y0
%where a=A, b=alpha, c=t, and y0=y0;
y=b(1)*exp(b(2)*x)+b(3);
end

