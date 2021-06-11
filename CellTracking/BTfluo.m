%BTfluo.m
%Rico Rojas, updated 1/21/19
%Calculates the average cytoplasmic fluorescence intensity from cell
%tracked with BacTrack.m.  

clear, close all

%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save fluorescent image stacks
%directories by themselves. 

%INPUT
%basename: name to save the results to.
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%icell: Cell array with length equal to the number of fluorescent
        %channels.  Each entry is a matrix (ncellxT) with the fluorescent intensities of each
        %cell, where rows are the cells and columns are time points.
%icell_av:  Cell array with length equal to the number of fluorescent
        %channels.  Each entry is a vector containing the population-
        %average of the single-cell fluorescent intensities.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='06062021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_aligned']; ['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_aligned']}; 
recrunch=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==0;

curdir=cd;
for i=1:length(channels)
    cd(channels{i}); 
    fluo_directory{i}=dir('*.tif');
end

cd(dirname)
load([basename '_BTphase'])
load([basename '_BTlab'])

intensities=cell(length(channels),1);


for i=1:length(channels)
    
    cd(channels{i}); 
    intensities_temp=zeros(size(lcell));
    
    for t=1:T
        t
        imagename=fluo_directory{i}(t).name;
        im=imread(imagename);
        for j=1:ncells
            intensities_temp(j,t)=mean(im(pixels{j,t}));
            
        end
    end
    intensities_temp(intensities_temp==0)=NaN;
    icell{i}=intensities_temp;
end

icell_av=cell(length(channels),1);
for i=1:length(channels)
    icell_av{i}=nanmean(icell{i});
end

save([basename '_BTfluo'])

elseif recrunch==1
    load ([basename '_BTfluo'])
end

%Plot data
%separate revelant time points
tidx=[1:2:T-1];

figure, hold on, 
plot(time(tidx),icell_av{1}(tidx), '-g')
plot(time(tidx),icell_av{2}(tidx), '-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
legend({'mNeonGreen', 'mCherry'})
saveas(gcf, [basename,'_intensity.fig'])
saveas(gcf, [basename,'_intensity.png'])


figure, hold on, 
for i=1:ncells
    plot(time(tidx),icell_av{1}(:,tidx), '-g')
end
%plot(time(tidx),icell_av{2}(tidx), '-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty

save([basename '_BTfluo'])