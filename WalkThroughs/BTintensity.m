%BTintensity.m  
%Calculates the intensity of a flourescent tracer dye that was imaged 
%simultaneously with a phase stack that was analyzed with BacTrack.m.

%INSTRUCTIONS FOR USE:
%Run BacTrack.m on a phase image stack.  Save the fluorescent image stack 
%(which images a tracer dye) in a separate directory.
%
%INPUT:
%basename: name of the image stack (from BacTrack).
%dirname: full path name of directory in which fluorescent image stack is
%           saved.
%Cmin: concentration of medium without tracer dye.
%Cmax: concentration of medium with tracer dye.
%
%OUTPUT:
%conc: vector of length T with concentration trace.

clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='032417_7_P0';%BS168 LB->LB3M->LB3M + 5%NLSS
basename='032417_7_P1';%BS168 LB->LB3M->LB3M + 5%NLSS
basename='032417_7_P2';%BS168 LB->LB3M->LB3M + 5%NLSS
basename='032417_7_P3';%BS168 LB->LB3M->LB3M + 5%NLSS
basename='032417_5_P0';%BS168 LB->LB3M
basename='032417_5_P1';%BS168 LB->LB3M
basename='032417_5_P2';%BS168 LB->LB3M
basename='032417_5_P3';%BS168 LB->LB3M
basename='032417_5_P4';%BS168 LB->LB3M
basename='032417_5_P5';%BS168 LB->LB3M
basename='040517_2_P0';%BS168 LB->LB+5% NLS
basename='040517_2_P1';%BS168 LB->LB+5% NLS
basename='040517_2_P2';%BS168 LB->LB+5% NLS
basename='040517_2_P3';%BS168 LB->LB+5% NLS
basename='040517_2_P4';%BS168 LB->LB+5% NLS
basename='071117_1_P0';%LB->LB+1M EDTA
basename='071117_1_P1';%LB->LB+1M EDTA
basename='071117_2_P0';%LB->LB+10M EDTA
basename='071117_2_P1';%LB->LB+10M EDTA
basename='071117_2_P2';%LB->LB+10M EDTA
basename='051915_3';%bacillus, membrane potential, 1M->LB
basename='080315_8_P0';skp=[];%LB500->LB400, SMB80
%basename='080315_8_P1';skp=[];%LB500->LB400, SMB80
%basename='080315_8_P2';skp=[];%LB500->LB400, SMB80
%basename='080315_8_P3';skp=[];%LB500->LB400, SMB80
%basename='080315_8_P4';skp=[];%LB500->LB400, SMB80
%basename='080315_8_P5';skp=[];%LB500->LB400, SMB80
basename='012820_1';skp=[];%fast 500 mM hypershock, Mg1655
basename='012920_016';skp=[];%fast 500 mM hypershock, Mg1655

Cmax=500;
Cmin=0;

recrunch=0;

%dirname=['/Volumes/FantomHD/Matlab Ready/' basename '/' basename '_fluo2_a'];
%dirname=['/Volumes/HD-PNTU3 1/Matlab Ready/' basename '/' basename '_fluo2_a'];
%dirname=['/Volumes/HD-PNTU3/Matlab Ready/' basename '/' basename '_fluo_a'];
%dirname=['/Volumes/HD-PNTU3/Matlab Ready/' basename '/' basename '_fluo'];
%dirname=['/Matlab Ready/' basename '/' basename '_fluo_a'];
%dirname=['/Volumes/Lacie/Matlab Ready/' basename '/' basename '_2_a'];
%dirname=['/Users/Rico/Documents/MATLAB/Matlab Ready/' basename '/' basename '_2_a'];
dirname=['/Users/Rico/Documents/MATLAB/Matlab Ready/' basename '/' basename '_3_a'];
%dirname=['/Users/Rico/Documents/MATLAB/Matlab Ready/' basename '/' basename '_fluo'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==1
    load([basename '_BTint'])
else

%load data from BacTrack.
load([basename '_BT'],'T','lsk','skp','time','tmid','vav','tpoints','lcell','time2','vstd','Leff')
path(dirname,path)

%time2=tpoints(end,1:T);

Tp=T-lsk;
int=zeros(1,Tp);

curdir=cd;
cd(dirname);
directory=dir('*.tif');
cd(curdir);

for t=1:Tp
    t
    %load image
    imagename=directory(t).name;
    im=imread(imagename);
    
    int(t)=mean(nonzeros(im));
end

%insert skipped frames
intf=int;
for i=1:lsk
     intf=[intf(1:skp(i)-1),(intf(skp(i))+intf(skp(i)-1))/2,intf(skp(i):end)];
end

%interpolate intensity at times at which cells were imaged
inti=interp1(time2,intf,time);
inti=intf;

%calculate concentration
Imax=max(inti);
Imin=min(inti);
a=(Cmax-Cmin)/(Imax-Imin);
b=Cmin-a*Imin;
conc=a*inti+b;

end

%plot results
figure
ciplot((vav-vstd)*3600,(vav+vstd)*3600,tmid,[0.75 0.75 1]),hold on
plot(tmid,vav*3600,'-r','LineWidth',1)
ylabel('Growth Rate (s^{-1})')
g = get(gca,'Position');
g(1:2) = g(1:2) + (1-0.9)/2*g(3:4);
g(3:4) = 0.9*g(3:4);
set(gca,'Position',g);
box off
apos=get(gca,'Position');
xtic=get(gca,'Xtick');
%ytic=get(gca,'Xtick');
lims=get(gca,'Xlim');
axes('Position',apos,'XTick',xtic,'Color','none','YAxisLocation','right',...
    'Xlim',lims);
hold on
plot(time2,conc)
xlabel('time (m)')
ylabel('C (mM)')
fig2pretty

figure,plot(conc(2:end),Leff)

%save
load([basename '_BT'])

% figure, title('Cell Length vs. Time')
% clf
% hold on
% for i=1:ncells
%     plot(time,lcell(i,:)) 
% end
% xlabel('time (s)')
% ylabel('length (/mum)')


%plot results
% figure
% plot(time(48:80),lcell(i,48:80),'-r','LineWidth',2)
% ylabel('Growth Rate (s^{-1})')
% g = get(gca,'Position');
% g(1:2) = g(1:2) + (1-0.9)/2*g(3:4);
% g(3:4) = 0.9*g(3:4);
% set(gca,'Position',g);
% box off
% apos=get(gca,'Position');
% xtic=get(gca,'Xtick');
% lims=get(gca,'Xlim');
% axes('Position',apos,'XTick',xtic,'Color','none','YAxisLocation','right',...
%     'Xlim',lims);
% hold on
% plot(time,conc)
% xlabel('time (s)')
% ylabel('C (mM)')

% figure
% plot(time(58:66),lcell(i,58:66),'-r','LineWidth',2)
% ylabel('Growth Rate (s^{-1})')
% g = get(gca,'Position');
% g(1:2) = g(1:2) + (1-0.9)/2*g(3:4);
% g(3:4) = 0.9*g(3:4);
% set(gca,'Position',g);
% box off
% apos=get(gca,'Position');
% xtic=get(gca,'Xtick');
% lims=get(gca,'Xlim');
% axes('Position',apos,'XTick',xtic,'Color','none','YAxisLocation','right',...
%     'Xlim',lims);
% hold on
% plot(time,conc)
% xlabel('time (s)')
% ylabel('C (mM)')

%fig2pretty

save([basename '_BTint'])
