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
basenames=["06022021_FITCK_e10_p005", "06022021_FITCK_e10_p008", "06022021_FITCK_e10_p010", "06022021_FITCK_e20_p005", "06022021_FITCK_e20_p008", "06022021_FITCK_e20_p010", "06022021_FITCK_e40_p005", "06022021_FITCK_e40_p008","06022021_FITCK_e40_p010"];
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/06022021_analysis/06022021_FITCK_figures'];
recrunch=0;
vis=0;
B=length(basenames);%number of main directories to analyze
exposure=[repelem(10,3), repelem(20, 3), repelem(40,3)];
psi=[5 8 10 5 8 10 5 8 10];
t1=[1:31];
t2=[174:355];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(dirname)
    load(['06202021_dm.mat'])
else
    
%go to directory where .mat files are stored
cd(dirname)

%pre-allocate variables
timescales=nan(B,1);
Ts=nan(B,1);
times1=cell(B,1);
times2=cell(B,1);
intensity1=cell(B,1);
intensity2=cell(B,1);
coeff1=nan(B,3);
coeff2=nan(B,3);
yhat1=cell(B,1);
yhat2=cell(B,1);
A=nan(B,1);
f=cell(B,1);
yhat2=cell(B,1);
thalf=nan(B,1);
mindiff=nan(B,1);
exp2=[];

%let's go through the mat files and store the variables we need
for b=1:B
    
    %initialize coeff
    coeff0=[40000, -0.1, 0.1];

    base=char(basenames(b))
    load([base '_pb'])

    timescales(b,1)=tscale;
    Ts(b,1)=T;
    times1{b,1}=time(t1);
    times2{b,1}=time(t2);
    intensity1{b,1}=intensityAvg(tp1);
    intensity2{b,1}=intensityAvg(tp2);
    A(b,1)=intensityAvg(1,1);
    
    %here, we can fit to a nonlinear exponential eqxn
    coeff1(b,:)=nlinfit(times1{b,1}, intensity1{b,1}, @exponential, coeff0);
    coeff2(b,:)=nlinfit(times2{b,1}, intensity2{b,1}, @exponential, coeff0);
    
    %calculate the predicted y values
    yhat1{b,1}=exponential(coeff1(b,:), times1{b,1});
    yhat2{b,1}=exponential(coeff2(b,:), times2{b,1});
    
    %plot to see how well the eqxn fits the data
    figure
    scatter(times1{b,1}, intensit1y{b,1}), hold on
    plot(times1{b,1}, yhat1{b,1})
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
%     
%     %calculate new yhat values
%     yhat2{b,1}=f{b,1}(times{b,1});
%     
    %pull out the coefficient values
    exp2(b,:)=coeffvalues(f{b,1});
end

save(['exposure_dm'])

end

%plot exposure as a function of intensity
k1=coeff(:,2);

figure
plot(exposure, k1)
ylabel('k1')
xlabel('exposure (ms)')
title('Coefficients as a function of Exposure')
saveas(gcf, [basename,'_exposureFunction.fig'])
saveas(gcf, [basename,'_exposureFunction.png'])

figure
plot(psi, k1)
ylabel('k1')
xlabel('exposure (ms)')
title('Coefficients as a function of Exposure')
saveas(gcf, [basename,'_psiFunction.fig'])
saveas(gcf, [basename,'_psiFunction.png'])

kb=exp2(:,2);
kd=exp2(:,4);

yyaxis left
plot(exposure, kb)
ylabel('b'), hold on
yyaxis right
plot(exposure,kd)
ylabel('d')
title('Coefficients as a function of Exposure')
saveas(gcf, [basename,'_exposureFunction2.fig'])
saveas(gcf, [basename,'_exposureFunction2.png'])

function [y] = exponential(b,x)
%this function calculates y=A*(e^alpha*t)+y0
%where a=A, b=alpha, c=t, and y0=y0;
y=b(1)*exp(-b(2)*x)+b(3);
end

