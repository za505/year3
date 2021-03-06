%photobleachMeasure.m
%Zarina Akbary, updated 08/02/21
%Calculates the rate of photobleaching

clear, close all

%INSTRUCTIONS FOR USE:
%save fluorescent image stacks directories by themselves

%INPUT
%basename: name of the file
%dirname: file location
%dye
%frameRate, in frames per second
%exposure, in ms
%flowRate, in psi
%flowTime, time interval during which dye is perfused in the chip, in
%minutes
%b=photobleaching constant

%OUTPUT:
%intensity=vector with raw intensity values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='07242021_Exp1';
photodir=['/Users/zarina/Documents/MATLAB/MatlabReady/07242021_analysis'];

dye=['FITCK'];
frameRate=120; %seconds but time array is in minutes
exposure=30;
flowRate=8;

recrunch=0;
vis=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savename)
    load([basename '_pb'])
else
    
cd(photodir);
directory=dir('*dm.mat');
N=length(directory);

%% Step 1: fit the data, find the photobleaching constant, and take the average
b=[];
cDiff=[];
count=0;
photoData=[];

for i=1:N
    
    load(directory(i).name)
    
    for n=1:height(norm_green)
        count=count+1;
        photoData(count,:)=norm_green(n,:);
        [xData, yData]=prepareCurveData(time, norm_green(n,:));
        f=fit(xData, yData, 'exp1'); %fit to a linear exponential
        if vis==1
            n
            plot(f, xData, yData), pause, close
        end
        
        b(count,1)=f.b; %extract the photobleaching constant
        cDiff(count,:)=diff(photoData(count,:)); %this is the dC*dt for the control
    end
end

bavg=abs(mean(b));

%% Step 2: integrate using the Newton Method 
dt=2;
bDiff=[];
dX=[];
X=zeros(height(photoData), length(photoData));
%dX/dt = dC/dt - BC
%X(n+1)=X(n) + dX/dt*dt
for i=1:height(photoData)
    bDiff(i,:) = bavg * photoData(i,2:end); 
    dX(i, :) = cDiff(i, :) + bDiff(i,:);
    X(i,1)=photoData(i,1);
    for n=1:length(photoData)-1
        X(i,n+1) = X(i,n) + dX(i,n)*dt;
    end
end

%smooth the data
for i=1:height(photoData)
    sX(i, :)=movingaverage(X(i, :),10);
end

%plot to check
if vis==1
     figure,hold on
    for i=1:height(photoData)
       
        plot(time, photoData(i,:))
        %plot(time, X(i,:))
%         plot(time, sX(i,:))
        %pause,close
    end
end

%% Step 3: perform this integration for the rest of the data
cd('../')
cd('mNeonGreenDiffusion_analysis/07202021_analysis') 

load('diffusionAnalysis.mat', 'fullData', 'time'); 

% %first, let's parse the data
cond1 = fullData(strcmp(fullData.Condition, "LB")==1, :);
cond2 = fullData(strcmp(fullData.Condition, "EDTA")==1, :);
cond3 = fullData(strcmp(fullData.Condition, "Mg2+")==1, :);

yData1=cond1(cond1.halfie==0,:);
yData2=cond2(cond2.halfie==0,:);
yData3=cond3(cond3.halfie==0,:);

experiments=unique(fullData.Experiment);

%first, LB
iData=yData1.("normalized intensity");
cDiff1=diff(iData, 1, 2);
bDiff1=[];
dX1=[];
X1=zeros(height(iData), length(iData));
for i=1:height(iData)
    bDiff1(i,:) = bavg * iData(i,2:end); 
    dX1(i, :) = cDiff1(i, :) - bDiff1(i,:);
        for n=1:length(iData)-1
            X1(i,n+1) = X1(i,n) + dX1(i,n)*dt;
        end
end

for i=1:height(yData1)
    figure,hold on
    plot(time, yData1.("normalized intensity")(i,:))
    plot(time, X1(i,:))
    pause,close
end

%EDTA
iData=yData2.("normalized intensity");
cDiff2=diff(iData);
bDiff2=[];
dX2=[];
X2=zeros(height(iData), length(iData));
for i=1:height(iData)
    bDiff2(i,:) = bavg * iData(i,2:end); 
    dX2(i, :) = cDiff2(i, :) - bDiff2(i,:);
        for n=1:length(iData)-1
            X2(i,n+1) = X2(i,n) + dX2(i,n)*dt;
        end
end
    
for i=1:height(yData1)
    figure,hold on
    plot(time, yData2.("normalized intensity")(i,:))
    plot(time, X2(i,:))
    pause,close
end

%Mg2+
iData=yData3.("normalized intensity");
cDiff3=diff(iData);
bDiff3=[];
dX3=[];
X3=zeros(height(iData), length(iData));
for i=1:height(iData)
    bDiff3(i,:) = bavg * iData(i,2:end); 
    dX3(i, :) = cDiff3(i, :) - bDiff3(i,:);
        for n=1:length(data)-1
            X3(i,n+1) = X3(i,n) + dX3(i,n)*dt;
        end
end

for i=1:height(yData1)
    figure,hold on
    plot(time, yData3.("normalized intensity")(i,:))
    plot(time, X3(i,:))
    pause,close
end
    

      
cd(savename)
save([basename '_pb'])

end

%plot figures
figure(1)
plot(time, intensityAvg)
title('Average Intensity vs Time')
xlabel('Time (s)')
ylabel('Average Intensity')
fig2pretty
xlim([-10 Inf])
saveas(gcf, [basename,'_intensityAvg.fig'])
saveas(gcf, [basename,'_intensityAvg.png'])

%%%%%Functions
% function [cDiff, dX, X]=NewtonMethod(norm_intensity, pbc) %pbc=photobleaching constant
%         
%     cDiff=diff(norm_intensity);
%     bDiff=[];
%     dX=[];
%     X=zeros(height(norm_intensity), length(norm_intensity));
%     for i=1:height(norm_intensity)
%         bDiff(i,:) = pbc * data(i,2:end); 
%         dX(i, :) = cDiff(i, :) - bDiff(i,:);
%             for n=1:length(data)-1
%                 X(i,n+1) = X(i,n) + dX(i,n)*dt;
%             end
%     end
% end 