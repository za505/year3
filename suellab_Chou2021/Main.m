clear all; %close all;


%parameters         
aQ = 0.01*0.001;    %Glutamine synthetase (GS) reaction rate (uM-2*h-1)
bQ = 1*0.001;       %GOGAT reaction rate (uM-2*h-1)
aF = 1;             %FBI-GS production rate through binding of glutamine metabolites to GS (uM*h-1)
bF = 1;             %FBI-GS degradation rate (h-1)
TT = 1;             %Total TnrA concentration (uM)
bT = 10;            %TnrA inhibition rate by FBI-GS (uM*h-1)
aA = 10;            %maximum GOGAT synthesis rate (uM*h-1)
bA = 1;             %GOGAT degradation rate (h-1)

n1 = 5;             %Hill coefficient for FBI-GS production
n2 = 5;             %Hill coefficient for the inhibition of GOGAT by TnrA

cq = 1;             %Glutamine degradation rate by metabolism (h-1)
aS = 10;            %Glutamate synthetase production rate (mM*h-1)
bS = 1;             %Glutamate synthetase degradation rate (h-1)
aE = 2*1000*10;     %Glutamate supply rate 10x (uM*h-1)
cE = 0.1;           %Glutamate degradation rate by other metabolic processes (h-1)
aN = 0.1;           %Ammonium production rate (h-1)
ao = 10*1000*10;    %2-oxoglutarate supply rate from the TCA cycle 10x (uM*h-1)
co = 10;            %2-oxoglutarate degradation rate through the TCA cycle (h-1)

KT = 0.1;           %Threshold of the inhibition of GOGAT by TnrA (uM)
KF = 1*1000;        %Threshold of FBI-GS production (uM)
KR = 0.3;           %Threshold transcriptional activation of PnasA by TnrA (uM)
aT = 1;             %TnrA activation rate by unbinding with FBI-GS (h-1)

aR = 10;            %maximum PnasA-yfp synthesis rate (uM*h-1)
bR = 1;             %YFP degradation rate (h-1)

r0 = 1.5;           %radius of area plated with initial B. subtitles culture
dr = 0.03;          %dr 
dt = 0.3;           %dt
eta = 2;            %?freeze factor? of metabolic rate
e = 0.13;           %metabolic scaling factor
cR = 0.05;          %YFP photobleaching rate (h-1)

k = dr/dt;          %biofilm expansion rate (mm*h-1)
%%
%Initial Conditions: [Q, F, T, A, S, E, N, O, R]
y0 = [0.6*1000; 0.7; 0.1; 2; 10; 90*1000; 1*1000; 10*1000; 0.1];

parameters = [aQ, bQ, aF, bF, TT, bT, aA, bA, n1, n2, cq, aS, bS, aE, aN, ao, co, KT, KF, aT, cE, aR, bR, KR, e];

% These are the cells at the edge of the biofilm   
[t,y] = ode45(@(t,y) NasA_osciIII_D(t,y,parameters),0:dt:dt*100,y0);
%%
% time-trace storage files Q_tc = zeros(position,time)
Q_tc = diag(y(2:end,1));
F_tc = diag(y(2:end,2));
T_tc = diag(y(2:end,3));
A_tc = diag(y(2:end,4));
S_tc = diag(y(2:end,5));
E_tc = diag(y(2:end,6));
N_tc = diag(y(2:end,7));
O_tc = diag(y(2:end,8));
R_tc = diag(y(2:end,9));
W_tc = zeros(100);


%for each position, simulate the entire timetrace
for x = 1:98
    
    t0 = 0;
    
    %Run for tfinal time. tfianl = (100 - x)*dt 
    tfinal = (100 - x)*dt;
    
    
    parameters = [aQ, bQ, aF, bF, TT, bT, aA, bA, n1, n2, cq, aS, bS, aE, aN, ao, co, KT, KF, aT, cE, aR, bR, KR, e, eta, cR];
    y0 = [Q_tc(x+1,x+1);F_tc(x+1,x+1);T_tc(x+1,x+1);A_tc(x+1,x+1);S_tc(x+1,x+1);E_tc(x+1,x+1);N_tc(x+1,x+1);O_tc(x+1,x+1);R_tc(x+1,x+1);0];
    
    % These are the cells left at position x 
    [t,y2] = ode45(@(t,y) NasA_osciIII_eta(t,y,parameters),0:dt:tfinal,y0);

    % Save the time trace at position x
    Q_tc(x,x+1:end) = y2(2:end,1);
    F_tc(x,x+1:end) = y2(2:end,2);
    T_tc(x,x+1:end) = y2(2:end,3);
    A_tc(x,x+1:end) = y2(2:end,4);
    S_tc(x,x+1:end) = y2(2:end,5);
    E_tc(x,x+1:end) = y2(2:end,6);
    N_tc(x,x+1:end) = y2(2:end,7);
    O_tc(x,x+1:end) = y2(2:end,8);
    R_tc(x,x+1:end) = y2(2:end,9);
    W_tc(x,x+1:end) = e./(1+y2(2:end,10));
end

Q_tc(Q_tc==0) = NaN;
F_tc(F_tc==0) = NaN;
T_tc(T_tc==0) = NaN;
A_tc(A_tc==0) = NaN;
S_tc(S_tc==0) = NaN;
E_tc(E_tc==0) = NaN;
N_tc(N_tc==0) = NaN;
O_tc(O_tc==0) = NaN;
R_tc(R_tc==0) = NaN;
W_tc(W_tc==0) = NaN;

%% Plotting the results in kymograph
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create surface
surface(R_tc(2:end,2:end),'Parent',axes1,'AlignVertexCenters','on','LineStyle', 'none');

ylabel({'position (a.u.)'});
xlabel({'time (a.u.)'});

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'Colormap',...
    [0 0 0.00520833348855376;0 0.0143222510814667 0.00498188426718116;0 0.0286445021629333 0.00475543504580855;0 0.0429667532444 0.00452898582443595;0 0.0572890043258667 0.00430253613740206;0 0.0716112554073334 0.00407608691602945;0 0.0859335064888 0.00384963769465685;0 0.100255757570267 0.00362318847328424;0 0.114578008651733 0.00339673925191164;0 0.1289002597332 0.00317029003053904;0 0.143222510814667 0.00294384057633579;0 0.157544761896133 0.00271739135496318;0 0.1718670129776 0.00249094213359058;0 0.186189264059067 0.00226449291221797;0 0.200511515140533 0.00203804345801473;0 0.214833766222 0.00181159423664212;0 0.229156017303467 0.00158514501526952;0 0.243478268384933 0.00135869567748159;0 0.2578005194664 0.00113224645610899;0 0.272122770547867 0.000905797118321061;0 0.286445021629333 0.000679347838740796;0 0.3007672727108 0.000452898559160531;0 0.315089523792267 0.000226449279580265;0 0.329411774873734 0;0 0.343137264251709 0;0 0.356862753629684 0;0 0.37058824300766 0;0 0.384313732385635 0;0 0.398039221763611 0;0 0.411764711141586 0;0 0.425490200519562 0;0 0.439215689897537 0;0 0.452941179275513 0;0 0.466666668653488 0;0 0.480392158031464 0;0 0.494117677211761 0;0 0.507843136787415 0;0 0.521568655967712 0;0 0.535294115543365 0;0 0.549019634723663 0;0 0.562745094299316 0;0 0.576470613479614 0;0 0.590196073055267 0;0 0.603921592235565 0;0 0.617647051811218 0;0 0.631372570991516 0;0 0.645098030567169 0;0 0.658823549747467 0;0 0.680147051811218 0;0 0.701470613479614 0;0 0.722794115543365 0;0 0.744117677211761 0;0 0.765441179275513 0;0 0.786764740943909 0;0 0.80808824300766 0;0 0.829411745071411 0;0 0.850735306739807 0;0 0.872058808803558 0;0 0.893382370471954 0;0 0.914705872535706 0;0 0.936029434204102 0;0 0.957352936267853 0;0 0.978676497936249 0;0 1 0],...
    'FontSize',12);
% Create colorbar
colorbar(axes1);
%%
filename = 'Test1_3xMSgg.mat';
save(filename);

%% Make movie for 1-D

outmovie2 = figure('Position',[300 300 600 400]);
set(gcf,'color','w');
VV = R_tc';
movie_name = '.mp4';
Movie = VideoWriter(movie_name, 'MPEG-4');
Movie.FrameRate = 20;  
open(Movie);


cmax = max(VV(end,:)*1.4);
%cmin = 0;
%cmax = max(max(UU))*1.4;
%cmax = 150;
cmin = 0;

for frame = 1:100
    plot(VV(frame,:), 'LineWidth', 2);%hold on;
    xlabel('Position (a.u.)','FontSize',12);
    ylabel('PnasA-yfp','FontSize',12);
%    xlim([0 100]);
%    set(gca, 'visible', 'off')
%    set(gca,'xtick',[]);
%    set(gca,'xticklabel',[]);  
%    set(gca,'ytick',[]);
%    set(gca,'yticklabel',[]);  
    ylim([cmin cmax]);
    xlim([0 100]);
    ax = gca;
%    set(ax,'YTick',100:100:500,'YTickLabel',100:100:500)
    ax.FontSize = 12;
    
    frames = getframe(outmovie2);
    writeVideo(Movie,frames);
end

close(Movie);
%% Plotting the kymograph into round biofilms
timepoint = 90;
R2 = repmat(R_tc(:, timepoint)', 61, 1);
figure
% Create polar data
[r,t] = meshgrid(r0+dr:dr:r0+100*dr,0:pi/120:(0.5*pi));
z = R2;
% Convert to Cartesian
xx = r.*cos(t);
yy = r.*sin(t);
h = polar(xx,yy);
hold on;
contourf(xx,yy,z,  100,'LineColor', 'none');
colormap([0 0 0.00520833348855376;0 0.0143222510814667 0.00498188426718116;0 0.0286445021629333 0.00475543504580855;0 0.0429667532444 0.00452898582443595;0 0.0572890043258667 0.00430253613740206;0 0.0716112554073334 0.00407608691602945;0 0.0859335064888 0.00384963769465685;0 0.100255757570267 0.00362318847328424;0 0.114578008651733 0.00339673925191164;0 0.1289002597332 0.00317029003053904;0 0.143222510814667 0.00294384057633579;0 0.157544761896133 0.00271739135496318;0 0.1718670129776 0.00249094213359058;0 0.186189264059067 0.00226449291221797;0 0.200511515140533 0.00203804345801473;0 0.214833766222 0.00181159423664212;0 0.229156017303467 0.00158514501526952;0 0.243478268384933 0.00135869567748159;0 0.2578005194664 0.00113224645610899;0 0.272122770547867 0.000905797118321061;0 0.286445021629333 0.000679347838740796;0 0.3007672727108 0.000452898559160531;0 0.315089523792267 0.000226449279580265;0 0.329411774873734 0;0 0.343137264251709 0;0 0.356862753629684 0;0 0.37058824300766 0;0 0.384313732385635 0;0 0.398039221763611 0;0 0.411764711141586 0;0 0.425490200519562 0;0 0.439215689897537 0;0 0.452941179275513 0;0 0.466666668653488 0;0 0.480392158031464 0;0 0.494117677211761 0;0 0.507843136787415 0;0 0.521568655967712 0;0 0.535294115543365 0;0 0.549019634723663 0;0 0.562745094299316 0;0 0.576470613479614 0;0 0.590196073055267 0;0 0.603921592235565 0;0 0.617647051811218 0;0 0.631372570991516 0;0 0.645098030567169 0;0 0.658823549747467 0;0 0.680147051811218 0;0 0.701470613479614 0;0 0.722794115543365 0;0 0.744117677211761 0;0 0.765441179275513 0;0 0.786764740943909 0;0 0.80808824300766 0;0 0.829411745071411 0;0 0.850735306739807 0;0 0.872058808803558 0;0 0.893382370471954 0;0 0.914705872535706 0;0 0.936029434204102 0;0 0.957352936267853 0;0 0.978676497936249 0;0 1 0])
% Create colorbar
%caxis([0 10]);
colorbar
% Hide the POLAR function data and leave annotations
set(h,'Visible','off')
% Turn off axes and set square aspect ratio
labels = findall(gca,'type','text');
set(labels,'visible','off');
labels = findall(gca,'type','line');
set(labels,'visible','off');
set(gcf,'color','w');
axis off
axis image
hold off;
