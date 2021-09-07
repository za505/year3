%LinTrack.m
%Tracks lineage of bacteria 
%To be used after BacTrack/2.m
%Dylan Fitzmaurice, last updated 033021
clear 
close all

%input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basename='07102021_Exp1';%Name of the image stack, used to save file.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07102021_analysis/' basename '_colony1/' basename '_phase2/' basename '_figures'];%Directory to save the output .mat file to.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(savedir)
load([basename '_BTphase'],'T','im','lcents','cellnum','tstamp','l','lscale', 'B', 'time', 'pixels', 'a');
load([basename '_BTlab'], 'labels');
ppxi=90;
%% Part 2: Lineage Tracking

%view binary movie
% for i=flip(1:T)
% imshow(labels(:,:,i))
% pause(0.1)
% clf
% end
%

%Get coordinates of cells
for i = flip(1:T)
numofcells(i)=max(max(labels(:,:,i)));% find number of cells
mxnumofcells=max(numofcells);
    for j=1:mxnumofcells
        [r1,c1] = find(bwlabel(labels(:,:,i))==j);   %Find coordinates of cells in frame T

        r1c1{j,i}=[r1 c1]; % concatenate XY coords frame T 
        
    end
end

for i=1:T
    for j=1:mxnumofcells
        if isempty(r1c1{j,i})==0 
            cell_id_na(j,i)=j; % cell id not aligned/ not lineage tracked
        end
    end
end

for i = flip(1:T)
    for j=1:mxnumofcells
        if isempty(r1c1{j,i})==1
            r1c1{j,i} = [NaN NaN];
        end
    end
end

%Align cells from frame to frame, correct r1c1 to track lineage
for i=T
    for j=1:mxnumofcells
        cell_id(j,i)=j; %makes cell id vector
        r2c2{j,i}=r1c1{j,i};%Seed r2c2
    end
end

%Tracks cells backward through time  
for i=flip(2:T)
    for j=1:mxnumofcells
        if isnan(r2c2{j,i})==1 %if NaN, find previous nonNaN id matrix for tracking
            index=cellfun('length', r2c2(j,:));
            nindex=index>1;
            ncols=min(find(nindex));
            Lia1=ismember(r2c2{j,ncols},r1c1{j,i-1},'rows');
            if 100*sum(Lia1)/length(r1c1{j,i-1})>ppxi==1
                cell_id(j,i-1)=j;
                r2c2{j,i-1}=r1c1{j,i-1};
            else
                cell_id(j,i-1)=NaN;
                r2c2{j,i-1}=NaN;
            end
        else
            nidx=cellfun('length', r2c2(j,:)); %Normal tracking backward 
            nnidx=nidx>1;
            ncols=min(find(nnidx));
            Lia1=ismember(r1c1{j,i-1},r2c2{j,ncols},'rows');
                if 100*sum(Lia1)/length(r2c2{j,ncols})>ppxi==1
                            cell_id(j,i-1)=j;
                            r2c2{j,i-1}=r1c1{j,i-1};
                    else
                        r2c2{j,i-1}=NaN;
                end
        end
     end
end

% Fill NaNs tracking backward
for i = flip(1:T)
for j=1:mxnumofcells
    for k=1:mxnumofcells
        if isnan(r2c2{j,i})==1
            index=cellfun('length', r2c2(j,:));
            nindex=index==1;
            ncols=max(find(nindex));
            Lia1=ismember(r2c2{j,ncols+1},r1c1{k,i},'rows');
            if 100*sum(Lia1)/length(r1c1{k,ncols+1})>ppxi==1
                cell_id(j,i)=k;
                r2c2{j,i}=r1c1{k,i};
            else
%                 cell_id(j,i-1)=NaN;
%                 r2c2{j,i-1}=NaN;
            end
        end
     end
end
end

%% Organize lcell with lineage tracking
%Track cells frame to frame
tracks=zeros(size(im));
[r2,c2] = find(bwlabel(labels(:,:,end))>=1);
linind=sub2ind(size(im),r2,c2);%new linind for tracking binary image
tracks(linind)=1;

[tracksL(:,:,i),ncells]=bwlabel(tracks);

nlcell=zeros(ncells,T);
nacell=zeros(ncells,T);
nB=cell(ncells,T);
npixels=cell(ncells,T);

%new lcell still original lineage tracking
for i=1:lcents
    nlcell(cellnum(i),tstamp(i))=l(cellnum(i),tstamp(i));
    nacell(cellnum(i),tstamp(i))=a(cellnum(i),tstamp(i));
    nB{cellnum(i),tstamp(i)}=B{cellnum(i),tstamp(i)};
    npixels{cellnum(i),tstamp(i)}=pixels{cellnum(i),tstamp(i)};
end

%
for i=1:T %length of cell id
    for j=1:mxnumofcells %number of cells
        for k=1:mxnumofcells %number of cells
            if cell_id(j,i)==cell_id_na(k,i)
                nnlcell(j,i)=nlcell(k,i);
                nnacell(j,i)=nacell(k,i);
                nnB{j,i}=nB{k,i};
                nnpixels{j,i}=npixels{k,i};
            end
        end
    end
end

%scale lcell
nnlcell=nnlcell*lscale;
nnacell=nnacell*lscale^2;

%%
%if cell lengths converge,divide points of same length in half
for i=1:T %length of cell id
    for j=1:height(nnlcell) %number of cells
        for k=1:height(nnlcell) %number of cells
            if j-k>=1 & nnlcell(j,i)==nnlcell(k,i)
                if nnlcell(j,i)==nnlcell(k,i)
                    nnlcell(j,i)=nnlcell(j,i)/2;
                    nnlcell(k,i)=nnlcell(k,i)/2;
                end
            end
        end
    end
end

tmid=(time(2:end)+time(1:end-1))/2;
%%
%plotting new lcell 
nnlcell(nnlcell == 0) = NaN;

% figure
plot(1:T,nnlcell)

%%Calculating the growth rate
deltat=time(2:end)-time(1:end-1);
nv=(nnlcell(:,2:end)-nnlcell(:,1:end-1))./((nnlcell(:,1:end-1)+nnlcell(:,2:end))/2);
av=(nnacell(:,2:end)-nnacell(:,1:end-1))./((nnacell(:,1:end-1)+nnacell(:,2:end))/2);
for i=1:height(nnlcell)
    nv(i,:)=nv(i,:)./deltat;
    av(i,:)=av(i,:)./deltat;
end

%%Generate plots
figure, title('Cell Length vs. Time')
clf
hold on
for i=1:height(nnlcell)  
    nnlcell(i,:)=movingaverage2(nnlcell(i,:),3);
    %indx=isnan(lcell(i,:))~=1;
    %indx=find(indx);
    %plot(time(indx),lcell(i,indx))
    plot(time(1:end),nnlcell(i,1:end)) 
end
xlabel('Time (s)')
ylabel('Length (\mum)')
fig2pretty
saveas(gcf,[basename,'_LinlTraces.png'])

figure, title('Elongation Rate vs. Time')
hold on
for i=1:height(nnlcell) 
    plot(tmid,nv(i,:))
end
%plot(tmid,vav,'-r')
xlabel('Time (s)')
ylabel('Elongation Rate (s^{-1})')
fig2pretty
saveas(gcf, [basename,'_LineTraces.png'])

%%Troubleshoot Check
for t=flip(T:1)
figure
imshow(im)
hold on
for k=1:height(nnB)
if isempty(nnB{k,t})==0
plot(nnB{k,t}(:,1), nnB{k,t}(:,2),'-r')
end
end
pause, close
end

%% Save data 
clear labels

cd(savedir)
save([basename '_LT'])

beep