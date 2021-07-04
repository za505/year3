%decayMeasure.m
%Zarina Akbary, updated 07/04/21
%Calculates changes in fluor. intensity. Incorporates BTfluo.m code.

clear, close all

%INSTRUCTIONS FOR USE:
%run BacTrack.m first

%INPUT
%basename: experiments of interest
%dirname: where .mat files are stored
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%timescale=array of time scales
%Ts=array of T (number of frames)
%switchF=array of switchFrames
%A=array of A values (A*e^(-alpha*t)+y0)
%f=cell of coeff for exponential eqxn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='06062021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_aligned']; ['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_aligned']}; 
recrunch=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(dirname)
    load([basename '_dm.mat'])
else
    
    for i=1:length(channels)
        cd(channels{i}); 
        fluo_directory{i}=dir('*.tif');
    end
    
    %go to directory where .mat files are stored
    cd(dirname)
    load([basename '_BTphase'], 'B', 'T', 'ncells', 'time', 'pixels')

    %pre-allocate variables
    icell_green=[];
    icell_mcherry=[];
    
    %calculate time index
    tidx=1:2:T-1; %remember, fluor images were taken every 2 frames
    
    for i=1:length(channels)
        
        cd(channels{i});
        intensity=nan(ncells, length(tidx));
        
        for j=1:length(tidx)
            t=tidx(j);
            
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            for n=1:ncells
                intensity(n,j)=mean(im(pixels{n,t}));
            end
            
        end
        
        for n=1:ncells
            if i==1 & mean(intensity(n, 1:5)>20000)
                icell_green=[icell_green; intensity(n,:)];
            elseif i==2 & mean(intensity(n, 1:5)>20000)
                icell_cherry=[icell_cherry; intensity(n,:)];
            else
                continue
            end
        end
    end
end

%% Plot data

%define ncells
ncells=height(boun);

%figure out which traces belongs to which channel
% figure, hold on
% for n=1:height(icell{1:1})
%     plot(time(tidx),icell{1:1}(n, tidx)) 
% end

nGreen=[1,3,6,7,15,17,18,22,23,24,25,26,32,33,35,38];
nCherry=[2,4,5,8,9,10,11,12,14,16,19,20,21,29,31,36,39];
%skipped 13,27,28,30,34, and 37 because over half the values are NaN

%split the cell intensity array into separate variables based on channel
icell_green=icell{1,1}(nGreen, tidx);
icell_mcherry=icell{1,2}(nCherry, tidx);

%identity time point where cells lyse
tpt=480;

%initialize coefficients
coeff0=[40000, 0.1];
    
for n=1:height(icell{1:1})
    n
    x=sum(isnan(icell{1,1}(n,tidx)))
end

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
    scatter(time, icell_green(i,:))
    %scatter(time, yhat1(i,:))
end
% saveas(gcf, [basename,'_expGreen.fig'])
% saveas(gcf, [basename,'_expGreen.png'])
 
figure, hold on
for i=1:height(icell_mcherry)
    plot(time, icell_mcherry(i,:))
    scatter(time, yhat2(i,:))
end
% saveas(gcf, [basename,'_expCherry.fig'])
% saveas(gcf, [basename,'_expCherry.png'])

%save(['dm'])

end
%%
% coeff3=nan(height(icell_mcherry), length(coeff0));
% yhat3=nan(height(icell_green), length(time));
% for i=1:height(icell_green)
%     coeff3(i,:)=nlinfit(time, icell_green(i,:), @linear, coeff0);
%     yhat3(i,:)=linear(coeff3(i,:), time);
% end
% 
% figure, hold on
% for i=1:height(icell_green)
%     plot(time, icell_green(i,:))
%     scatter(time, yhat3(i,:))
% end
% saveas(gcf, [basename,'_expGreen2.fig'])
% saveas(gcf, [basename,'_expGreen2.png'])

%plot to see how well the eqxn fits the data
figure(1), hold on
for i=1:height(icell_green)
    plot(time(tidx), icell_green(i,:))
    %scatter(time, yhat1(i,:))
end
xline(tpt, '--', {'Membrane Lysis'})
title('Cellular Intensity of mNeonGreen vs Time')
xlabel('Time (s)')
ylabel('Cellular Intensity (A.U.)')
saveas(gcf, [basename,'_fullGreen.fig'])
saveas(gcf, [basename,'_fullGreen.png'])
 
figure(2), hold on
for i=1:height(icell_mcherry)
   plot(time(tidx), icell_mcherry(i,:))
end
xline(tpt, '--', {'Membrane Lysis'})
title('Cellular Intensity of mCherry vs Time')
xlabel('Time (s)')
ylabel('Cellular Intensity (A.U.)')
saveas(gcf, [basename,'_fullCherry.fig'])
saveas(gcf, [basename,'_fullCherry.png'])

intensity=cell(ncells,T-1);

for n=1:ncells
    
    if ismember(n, nGreen)
        [~, locB]=ismember(n, nGreen);
        cd( strcat(dirname, '/', basenames(1), '/', basename, '_figures'))
        
        dx=max(boun{n,1}(:,1))-min(boun{n,1}(:,1)); %columns are the x direction
        dy=max(boun{n,1}(:,2))-min(boun{n,1}(:,2)); %rows are the y direction
        x=min(boun{n,1}(:,1))-dx/2:max(boun{n,1}(:,1))+dx/2;
        y=min(boun{n,1}(:,2))-dy/2:max(boun{n,1}(:,2))+dy/2;
        [X,Y] = meshgrid(x,y);
        
        for t=1:2:T-1
            t
            intensity{n,t}=ones(size(boun{n,t},1),2);
            intensity{n,t}=intensity{1,t}*65530;
        end

        v = VideoWriter('pt2_cyto','MPEG-4');
        open(v);

        figure, hold on
        for t=1:2:T
    plot3(boun{1,t}(:,1),boun{1,t}(:,2),intensity2{1,t})
    hold on
    t=surf(X,Y,cyto{1,t})
    rotate3d on;
    t.EdgeColor = 'interp';
    t.FaceColor = 'interp';
    view(0,95)
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause
    clf
end

close(v)
close all

    elseif ismember(n, nCherry)
        [~, locB]=ismember(n, nCherry);
        cd( strcat(dirname, '/', basenames(2), '/', basename, '_figures'))
    else
        continue
    end
end
        
cd(savedir)
% % v = VideoWriter('pt2_cyto','MPEG-4');
% % open(v);
% 

% 
% figure
% hold on
% for t=1:T
%     plot3(boun{1,t}(:,1),boun{1,t}(:,2),intensity2{1,t})
%     hold on
%     t=surf(X,Y,cyto{1,t})
%     rotate3d on;
%     t.EdgeColor = 'interp';
%     t.FaceColor = 'interp';
%     %view(0,95)
% %     frame = getframe(gcf);
% %     writeVideo(v,frame);
%     pause
%     clf
% end
% 
% %close(v)
close all

%%%%%%%%%%%Functions
function [y] = exponential(b,x)
%this function calculates y=A*(e^-t/tau)
%where b(1)=A, b(2)=tau, x=t, and the cellular intensity=y;
y=b(1)*exp(-x./b(2));
end

function [y] = linear(b,x)
%this function calculates y=A*(e^-t/tau)
%where b(1)=A, b(2)=tau, x=t, and the cellular intensity=y;
y=x.*b(1)+b(2);
end