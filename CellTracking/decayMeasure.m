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
t2=[174:357];
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

%initialize coeff
coeff0=[40000, -0.1, 0.1];
    
%let's go through the mat files and store the variables we need
for b=1:B

    base=char(basenames(b))
    load([base '_pb'])

    timescales(b,1)=tscale;
    Ts(b,1)=T;
    times1{b,1}=time(t1);
    times2{b,1}=time(t2)-(t2(1)-1);
    intensity1{b,1}=intensityAvg(t1);
    intensity2{b,1}=intensityAvg(t2);
    A(b,1)=intensityAvg(1,1);
    
    %here, we can fit to a nonlinear exponential eqxn
    coeff1(b,:)=nlinfit(times1{b,1}, intensity1{b,1}, @exponential, coeff0);
    coeff2(b,:)=nlinfit(times2{b,1}, intensity2{b,1}, @exponential, coeff0);
    
    %calculate the predicted y values
    yhat1{b,1}=exponential(coeff1(b,:), times1{b,1});
    yhat2{b,1}=exponential(coeff2(b,:), times2{b,1});
    
    %plot to see how well the eqxn fits the data
    figure
    scatter(times1{b,1}, intensity1{b,1}), hold on
    plot(times1{b,1}, yhat1{b,1}), title("Fit for k1")
    pause
    saveas(gcf, [basename,'_exp1.fig'])
    saveas(gcf, [basename,'_exp1.png'])
    close
    
    %plot to see how well the eqxn fits the data
    figure
    scatter(times2{b,1}, intensity2{b,1}), hold on
    plot(times2{b,1}, yhat2{b,1}), title("Fit for k2")
    pause
    saveas(gcf, [basename,'_exp2.fig'])
    saveas(gcf, [basename,'_exp2.png'])
    close
end

save(['dm'])

end

%plot exposure as a function of intensity
k1=coeff1(:,2);
k2=coeff2(:,2);

figure, hold on
plot(exposure([1 4 7]), k1([1 4 7]), 'r')
plot(exposure([2 5 8]), k1([2 5 8]), 'b')
plot(exposure([3 6 9]), k1([3 6 9]), 'g')
ylabel('k1')
xlabel('exposure (ms)')
title('Coefficients as a function of Exposure')
legend({'5 psi', '8 psi', '10 psi'})
saveas(gcf, [basename,'_exposureFunction.fig'])
saveas(gcf, [basename,'_exposureFunction.png'])

figure, hold on
plot(psi(1:3), k1(1:3), 'r')
plot(psi(4:6), k1(4:6), 'b')
plot(psi(7:9), k1(7:9), 'g')
ylabel('k1')
xlabel('flow rate (psi)')
legend({'10 ms', '20 ms', '40 ms'})
title('Coefficients as a function of Flow Rate')
saveas(gcf, [basename,'_psiFunction.fig'])
saveas(gcf, [basename,'_psiFunction.png'])

figure, hold on
plot(exposure([1 4 7]), k2([1 4 7]), 'r')
plot(exposure([2 5 8]), k2([2 5 8]), 'b')
plot(exposure([3 6 9]), k2([3 6 9]), 'g')
ylabel('k2')
xlabel('exposure (ms)')
title('Coefficients as a function of Exposure')
legend({'5 psi', '8 psi', '10 psi'})
saveas(gcf, [basename,'_exposureFunction2.fig'])
saveas(gcf, [basename,'_exposureFunction2.png'])

figure, hold on
plot(psi(1:3), k2(1:3), 'r')
plot(psi(4:6), k2(4:6), 'b')
plot(psi(7:9), k2(7:9), 'g')
ylabel('k2')
xlabel('flow rate (psi)')
legend({'10 ms', '20 ms', '40 ms'})
title('Coefficients as a function of Flow Rate')
saveas(gcf, [basename,'_psiFunction2.fig'])
saveas(gcf, [basename,'_psiFunction2.png'])

function [y] = exponential(b,x)
%this function calculates y=A*(e^alpha*t)+y0
%where a=A, b=alpha, c=t, and y0=y0;
y=b(1)*exp(b(2).*x)+b(3);
end

