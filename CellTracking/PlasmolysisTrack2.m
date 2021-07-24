%PlasmolysisTrack2.m
%Tracks bacterial growth from phase image stacks.  
%Customized for B. subtilis.
%Bactrack2.m is customized for B. subtilis filaments.


%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save phase image stack in a directory
%by itself.  Also save the micromanager meGFPta file as 'basename.txt' in
%the matlab path.
%
%INPUT:
%basename: name of the image stack.
%dirname:the full pathname of the directory where you saved the image
%        stack.
%metaname(optional):full or relative pathname of micromanager meGFPta file from
%which to extract time points.  If it is relative path name, the
%directory in which it is saved must be on the matlab path.
%lscale: microscope calibration in microns per pixels.
%sm: width of the Gaussian filter used in edge finder equals sm*sqrt(2).
%minL: minimum length of cells;
%minW: minimum width of cells;
%maxW: maximum width of cells;
%recrunch:0 or 1.  if you've already tracked the data set and just want to
%         re-plot the data enter 1.
%
%OUTPUT:
%T: number of time points.
%time: vector of length T with time points.
%tmid: vector of length T-1 with interstitial time points.
%ncells: number of individual cells tracked.
%lcell: ncells x T matrix of cell lengths.
%wcell: ncells x T matrix of cell widths.
%acell: ncells x T matrix of cell areas
%ew: ncells x T matrix of circumferential strains.
%acell: ncells x T matrix of cell areas.
%v: ncells x T-1 matrix of cell strain rates.
%B: ncells x T cell array with cell contours.
%mlines: ncells x T cell array with cell midlines
%wav: vector of length T with average cell widths.
%wstd: vector of length T with standard deviations of cell widths.
%wste: vector of length T with standard error of cell widths.
%vav: vector of length T-1 with average cell strain rate.
%vstd: vector of length T-1 with standard deviation of strain rates.
%vste: vector of length T-1 with standard error of strain rates.
%avav: vector of length T-1 with average cell areal strain rate.
%avstd: vector of length T-1 with standard deviation of areal strain rates.
%avste: vector of length T-1 with standard error of areal strain rates.
%ndp: vecotr of lenth T-1 with number of data points averaged over.

%Calls on the following m-files:
%norm16bit.m
%polefinder.m
%cellcurvature.m
%meGFPta.m
%extrema.m
%EffectiveLength.m
%fig2pretty.m
%movingaverage.m

clear
close all

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='05262021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05262021_analysis/' basename '/05262021_plasmolysisTrack2/05262021_TADA'];%Directory that the image stack is saved in.
phasename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05262021_analysis/' basename '/05262021_plasmolysisTrack2/05262021_phase'];
cytoname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05262021_analysis/' basename '/05262021_plasmolysisTrack2/05262021_GFP'];
savedir=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05262021_analysis/' basename '/05262021_plasmolysisTrack2/05262021_figures'];%Directory to save the output .mat file to.
%metaname=['/Users/Rico/Documents/MATLAB/Matlab Ready/' basename '/meGFPta.txt'];%Name of meGFPta file.  Will only work if images were taken with micromanager.
lscale=0.08;%%Microns per pixel.
multiScale=0;
tscale=10;%Frame rate.
% tscale2=1;
% tpt1=120; %number of seconds passed by first time set
% tpt2=240; %number of seconds passed by second time set
% tpt3=480; %number of seconds passed by third time set
% tpt4=1320; %number of seconds passed by fourth time step
thresh=0;%For default, enter zero.
IntThresh=20000;%Threshold used to enhance contrast. Default:35000
dr=1;%Radius of dilation before watershed 
sm=2;%Parameter used in edge detection %default sm=2
minL=2;%Minimum cell length
minW=0.2;%Minimum cell width
maxW=1.5;%Maximum cell width
minA=50;%Minimum cell area. default 50
cellLink=4;%Number of frames to ignore missing cells when tracking frame to frame
recrunch=0;%Display data from previously crunched data? 0=No, 1=Yes.
vis=1;%Display cell tracking? 0=No, 1=Yes.
checkhist=0;%Display image histogram? 0=No, 1=Yes.
hT=1; %hard code T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% let's load the phase images (note that the frames in phase should correspond to the frames in TADA)
cd(savedir)
load([basename '_WT'])

%Determine number of frames
cd(cytoname);
cytodir=dir('*.tif');
imM=imm;
imN=imn;

cim=nan(imM,imN);

percent_plas_cyto=cell(ncells,T);
plasm_cyto=nan(ncells, T);
percent_plasmolysis_cyto=nan(ncells,T);
cyto=cell(1,T);

for t=1:T
    t
    
    close all
    cd(cytoname);
    
    %Load phase image
    imagename=cytodir(t).name;
    img=imread(imagename);
    cyto{1,t}=img;

    %Normalize images
    ppix=0.5;
    img=norm16bit(img,ppix);
    
    if isempty(B{n,t})==0
    %overlay TADA boundaries
    figure
    imshow(img, [])
    hold on
    for n=1:ncells
        if isempty(B{n,t})==0
            plot(B{n,t}(:,1),B{n,t}(:,2),'-r')
        end
    end
    pause, close
    
    %threshold so lower pixel values are set to 0
    gfpThresh=42000;
    
    for j=1:imN
        for i=1:imM
            if img(i,j)<gfpThresh
                cim(i,j)=0;
            elseif img(i,j)>gfpThresh
                cim(i,j)=img(i,j);
            else
                continue
            end
        end
    end
  
    cim=uint16(cim);
    
    figure
    imshow(cim, [])
    hold on
    for n=1:ncells
        if isempty(B{n,t})==0
            plot(B{n,t}(:,1),B{n,t}(:,2),'-r')
        end
    end
    pause, close
    %cd(savedir), saveas(gcf, [basename,'_green_' num2str(t) '.png']), close
   
    imshow([img, cim])
    %saveas(gcf, [basename, '_img_cim_' num2str(t) '.png'])
    pause, close
    
    for n=1:ncells
        
        percent_plas_cyto{n,t}=double(cim)-A_in{n,t};
        
        for j=1:imN %x direction
            for i=1:imM %y direction
                if A_in{n,t}(i,j)==0
                    percent_plas_cyto{n,t}(i,j)=NaN;
                end
            end
        end
        
        for j=1:imN %x direction
            for i=1:imM %y direction
                if percent_plas_cyto{n,t}(i,j)<0
                    percent_plas_cyto{n,t}(i,j)=NaN;
                elseif percent_plas_cyto{n,t}(i,j)>0
                    percent_plas_cyto{n,t}(i,j)=1;   
                else
                    continue
                end
            end
        end
        
        plasm_cyto(n,t)=sum(nansum(percent_plas_cyto{n,t}(:,:)));
        percent_plasmolysis_cyto(n,t)=100-(plasm_cyto(n,t)/total_pix(n,t))*100;
    end

    for n=1:ncells
        for j=1:imN %x direction
            for i=1:imM %y direction
                if A_in{n,t}(i,j)==0
                    cyto{1,t}(i,j)=0;
                end
            end
        end
    end
    
    imshow(cyto{1,t}, []), pause
    else
        continue
    end
end

cd(savedir)
% % v = VideoWriter('pt2_cyto','MPEG-4');
% % open(v);
% 
% for t=1:T
%     intensity2{1,t}=ones(size(boun{1,t},1),2);
%     intensity2{1,t}=intensity{1,t}*2000;
% end
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

cd(savedir)
%save([basename '_PTgreen'])