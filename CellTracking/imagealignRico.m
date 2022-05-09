%ImageAlign.m
%Registers image stacks using a common reference stack.  
%Calls on dftregistration.m
%
%INSTRUCTIONS FOR USE:
%Save each image stack you want to register in separate directories.  Also
%save the reference image stack (which contains features from which it is 
%possible to track drift) in its own directory.  The program will create
%new directories for each of the image stacks to be registered in the same
%parent directory.  The reference stacks and the stacks to be registered
%should have the same number of images.
%
%INPUT
%basename: choose a name for the directories in which the aligned images
%          will be saved.
%dirname:  cell array containing the paths of the directories in which the
%          image stacks to be registered are.
%regname:  path of reference directory.
%
%CALLS ON:
%dftregistration.m
%imtranslate.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for 05062022_Exp1_colony1
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='05062022_Exp1';
dirname={['/Users/zarina/Downloads/NYU/Year3_2022_Spring/05062022_analysis/' basename '_phase/' basename '_full']; ['/Users/zarina/Downloads/NYU/Year3_2022_Spring/05062022_analysis/' basename '_mNeonGreen/' basename '_full']}; %['/Users/zarina/Downloads/NYU/Year3_2022_Spring/04052022_analysis/' basename '/'   basename '_colony4/' basename '_TADA/' basename '_full'];};
regname=['/Users/zarina/Downloads/NYU/Year3_2022_Spring/05062022_analysis/' basename '_phase/' basename '_crop']; 
%frame=[11,18]; %plasmolyzed frame that is aligned
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

workdir=pwd;
ld=length(dirname); %number of channels being aligned

cd(regname)
directory=dir('*.tif');
path(regname,path)
T=length(directory); %number of frames in each image stack
tf=T*ld; %number of total frames to be aligned

imagename=directory(1).name; %load the first frame in the register image stack
RefIm=imread([regname '/' imagename]);

count=0; %to keep track of how many frames have been aligned
for i=1:ld
    
    cd(dirname{i})
    path(dirname{i},path)
    directory2=dir('*.tif');
    cd('../')
    
    mkdir([basename '_aligned'])
    cd(['./' basename '_aligned'])
    
    shft=zeros(T,2);
    
    for t=1:T
        count=count+1;
        pd=count/tf*100;
        pds=sprintf(['%4.1f'],pd);
        [pds '%'] 

        imagename=directory(t).name;
        TestIm=imread([regname '/' imagename]);

        [output NewImFT]=dftregistration(fft2(RefIm),fft2(TestIm),10); %all subsequent images are aligned to the first one in the stack
        NewIm=abs(ifft2(NewImFT));

        shft(t+1,1)=output(3);
        shft(t+1,2)=output(4);

        imagename=directory2(t).name;

        I=imread([dirname{i} '/' imagename]);
        [counts,bins]=imhist(I);
        [~,maxpos]=max(counts);
        padcolor=bins(maxpos);

        I=imtranslate(I,shft(t+1,:),padcolor);

        b=sprintf(['%4.4d'],t);  
        savename=[basename '_a' b '.tif'];

        imwrite(I,savename);
    end
end



