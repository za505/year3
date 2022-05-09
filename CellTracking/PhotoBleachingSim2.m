%% Modeling Photobleach Correction
%Authors: Enrique Rojas and Zarina Akbary
%Date: 04/06/2022
%Purpose: to model the photobleach correction

clear, close all

%% User input
colorcode = {'#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#000000'};
savedir = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/04072022_analysis';
cd(savedir)

%% Simulate Diffusion and Photobleaching (at a 1-minute frame rate)
%define the time vector
N=240; %the number of minutes in 4 hours
dt=1;
time=[0:N]*dt;

%simulating the loss of fluorescence due solely to diffusion of
%fluorophore,which is dependent on the permeability time constant tau. A is
%the initial fluorescence value (normalized to 1), B is the final
%fluorescence value (set to 0, for simplicity), and t is time.
diffusion = @(t, A, B, tau)(A * exp(-t./tau) + B);
diffusionDerv = @(t, A, B, tau)(-A/tau * exp(-t./tau));

tau=[10, 50, 100, 200, 1000];
labels = {'10', '50', '100', '200', '1000'};

%pre-allocate variables
nrow = length(tau);
C = nan(nrow, N+1);

%simulate diffusion
for i=1:nrow
    C(i, :) = diffusion(time, 1, 0, tau(i));
end

%plot this diffusion trace over time
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',20), hold on
%subplot(1, 2, 1), hold on
for i=1:nrow
    plot(time, C(i, :), 'Color', colorcode{i}, 'LineWidth', 1.5)
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%ylim([0 Inf])
lgd = legend(labels)
title(lgd,'\tau')
title('Diffusion Simulation')
%saveas(gcf, 'simDiffusion.png')

% 'Add' Photobleaching the the Fluorescence Trace
%pre-allocate a vector to store the value of unbleached fluor. (Cu)
Cu = nan(nrow, N+1);
dCu = nan(nrow,N);
Cu(:,1)=1; %initialize the vector to 1, the normalized fluor. value
fr = dt; %the frame rate
beta = 33.5358; %the slope used to calculate alpha
b =  -0.0204; %the y-intercept

%now 'add' photobleaching to the trace, using the Newton Method
photobleaching = @(Cu, alpha) Cu/alpha;
alpha = (beta * fr) + b; 

for i=1:nrow
    for j=1:N
        dCb = photobleaching(Cu(i, j), alpha);
        dCp = diffusionDerv(time(j), 1, 0, tau(i));
        dCu(i, j) = (-dCb + dCp) * dt;
        Cu(i, j+1) = Cu(i, j) + dCu(i, j);
    end
end

%the wrong way to simulate photobleaching
% Cw=C;
% for i=1:N
%     dCb = photobleaching(C(i), alpha);
%     Cw(i+1)=C(i)-dCb;
% end

%plot the 'true' and photobleached fluor. traces
subplot(1, 2, 2), hold on
for i=1:nrow
    plot(time, C(i, :), 'Color', colorcode{i}, 'LineWidth', 1.5)
    plot(time, Cu(i, :), '--', 'Color', colorcode{i}, 'LineWidth', 1.5)
end
%plot(time, Cw, '--', 'Color', colorcode{8}, 'LineWidth', 1.5)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
title('Photobleaching Simulation (1-minute Frame Rate)')
%ylim([0 Inf])
%saveas(gcf, 'simPhotobleaching.png')

% 'Add' Photobleaching the the Fluorescence Trace at a 5-minute frame rate
%pre-allocate a vector to store the value of unbleached fluor. (Cu)
Cu = nan(nrow, N+1);
dCu = nan(nrow,N);
Cu(:,1)=1; %initialize the vector to 1, the normalized fluor. value
fr = 5; %the frame rate
beta = 33.5358; %the slope used to calculate alpha
b =  -0.0204; %the y-intercept

%now 'add' photobleaching to the trace, using the Newton Method
photobleaching = @(Cu, alpha) Cu/alpha;
alpha = (beta * fr) + b; 

for i=1:nrow
    for j=1:N
        dCb = photobleaching(Cu(i, j), alpha);
        dCp = diffusionDerv(time(j), 1, 0, tau(i));
        dCu(i, j) = (-dCb + dCp) * dt;
        Cu(i, j+1) = Cu(i, j) + dCu(i, j);
    end
end

%plot this diffusion trace over time
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',20)
subplot(1, 2, 1), hold on
for i=1:nrow
    plot(time, C(i, :), 'Color', colorcode{i}, 'LineWidth', 1.5)
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%ylim([0 Inf])
lgd = legend(labels)
title(lgd,'\tau')
title('Diffusion Simulation')

%plot the 'true' and photobleached fluor. traces
subplot(1, 2, 2), hold on
for i=1:nrow
    plot(time, C(i, :), 'Color', colorcode{i}, 'LineWidth', 1.5)
    plot(time, Cu(i, :), '--', 'Color', colorcode{i}, 'LineWidth', 1.5)
end
%plot(time, Cw, '--', 'Color', colorcode{8}, 'LineWidth', 1.5)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
title('Photobleaching Simulation (5-minute Frame Rate)')
%ylim([0 Inf])
%saveas(gcf, 'simPhotobleaching.png')
%% Regenerate the Original Trace Using the Current Correction Method
photobleached = @(Cu, dt, beta, b)(Cu/(beta*dt+b))*dt; %what portion of the unbleached trace is getting bleached?

%pre-allocate vectors
Cnew = nan(1, N+1); %the corrected trace
unbFrac = nan(1, N+1); %the fractiono of unbleached fluor.
Cbl = nan(1, N+1); %the amount of bleached fluor.

Cnew(1) = Cu(1);
unbFrac(1) = 1;
Cbl(1) = 0;
   
for i=1:N

    %dt = time(i+1) - time(i);
    dCB = photobleached(Cu(i), dt, beta, b); %this is the amount of photobleaching that occured in the bleached trace
    dCT = Cu(i+1) - Cu(i); %this is the total fluor. loss for the bleached trace
    dCP = dCT + dCB; %this is the amount of loss attributable to permeability
    dCP = dCP * unbFrac(i); %Correcting for the fact that a fraction of fluorophores are unbleached
    Cnew(i+1) = Cnew(i) + dCP;
    Cbl(i+1) = Cbl(i) + dCB + dCP * (1-unbFrac(i)); %Accounting fo the change in concentration fo bleached fluorophores
    unbFrac(i+1)=(Cu(i+1))/(Cu(i+1)+Cbl(i+1)); %Calculate the new fraction of unbleached fluorophores

end    

figure, hold on
plot(time, C, 'Color', colorcode{5}, 'LineWidth', 1.5)
plot(time, Cu, 'Color', colorcode{2}, 'LineWidth', 1.5)
plot(time, Cnew, 'LineStyle', '--', 'Color', colorcode{7}, 'LineWidth', 1.5)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
ylim([0 Inf])
saveas(gcf, 'simCurrentCorrection.png')

%% Regenerate the Original Trace Using the Previous Correction Method
%pre-allocate vectors
Cnew = nan(1, N+1); %the corrected trace
unbFrac = nan(1, N+1); %the fractiono of unbleached fluor.
Cbl = nan(1, N+1); %the amount of bleached fluor.

Cnew(1) = Cu(1);
unbFrac(1) = 1;
Cbl(1) = 0;

for i=1:N

    %dt = time(i+1) - time(i);
    dCB = photobleached(Cu(i), dt, beta, 0); %this is the amount of photobleaching that occured in the bleached trace
    dCT = Cu(i+1) - Cu(i); %this is the total fluor. loss for the bleached trace
    dCP = dCT + dCB; %this is the amount of loss attributable to permeability
    Cnew(i+1) = Cnew(i) + dCP;
    Cbl(i+1) = Cbl(i) + dCB + dCP * (1-unbFrac(i)); %Accounting fo the change in concentration fo bleached fluorophores
    unbFrac(i+1)=(Cu(i+1))/(Cu(i+1)+Cbl(i+1)); %Calculate the new fraction of unbleached fluorophores

end 

figure, hold on
plot(time, C, 'Color', colorcode{5}, 'LineWidth', 1.5)
plot(time, Cu, 'Color', colorcode{2}, 'LineWidth', 1.5)
plot(time, Cnew, 'LineStyle', '--', 'Color', colorcode{7}, 'LineWidth', 1.5)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
ylim([0 Inf])
saveas(gcf, 'simOriginalCorrection.png')

%% Try Modifying the Current Method to Regenerate the Trace
photobleached = @(Cu, dt, beta, b)(Cu/(beta*dt+b))*dt; %what portion of the unbleached trace is getting bleached?

%pre-allocate vectors
Cnew = nan(1, N+1); %the corrected trace
unbFrac = nan(1, N+1); %the fractiono of unbleached fluor.
Cbl = nan(1, N+1); %the amount of bleached fluor.

Cnew(1) = Cu(1);
unbFrac(1) = 1;
Cbl(1) = 0;
   
for i=1:N
    
    if unbFrac(i)==1
        dCB = photobleached(Cu(i), dt, beta, b); %this is the amount of photobleaching that occured in the bleached trace
        dCT = Cu(i+1) - Cu(i); %this is the total fluor. loss for the bleached trace
        dCP = dCT + dCB; %this is the amount of loss attributable to permeability
        Cbl(i+1) = Cbl(i) + dCB; %Accounting fo the change in concentration fo bleached fluorophores
        Cnew(i+1) = Cnew(i) + dCP;
        unbFrac(i+1)=(Cu(i+1))/(Cu(i+1)+Cbl(i+1)); %Calculate the new fraction of unbleached fluorophores
    elseif unbFrac(i)==0
        dCB = photobleached(Cu(i), dt, beta, b); 
        dCT = Cu(i+1) - Cu(i); 
        dCP = dCT + dCB; 
        Cbl(i+1) = Cbl(i) + dCB + dCP; 
        Cnew(i+1) = Cnew(i) + dCP;
        unbFrac(i+1)=(Cu(i+1))/(Cu(i+1)+Cbl(i+1)); 
    else
        dCB = photobleached(Cu(i), dt, beta, b); 
        dCT = Cu(i+1) - Cu(i); 
        dCP = dCT + dCB; 
        Cbl(i+1) = Cbl(i) + dCB + (dCP / (1-unbFrac(i))); 
        dCP = dCP / unbFrac(i); %Correcting for the fact that a fraction of fluorophores are unbleached
        Cnew(i+1) = Cnew(i) + dCP;
        unbFrac(i+1)=(Cu(i+1))/(Cu(i+1)+Cbl(i+1));         
    end
end    

figure, hold on
plot(time, C, 'Color', colorcode{5}, 'LineWidth', 1.5)
plot(time, Cu, 'Color', colorcode{7}, 'LineWidth', 1.5)
plot(time, Cnew, 'LineStyle', '--', 'Color', colorcode{8}, 'LineWidth', 1.5)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
ylim([0 Inf])
pause 

%% Now, start with a photobleached trace and then add permeability. Does the correction method still work? 
%define the time vector
N=3600;
dt=0.0167;
time=[0:N]*dt;
lysis = 7; %when diffusion starts
idx = sum(time<=lysis);

%define new 'true' trace
C = nan(1, N+1);
C(time<=lysis) = 1;
C(time>lysis) = diffusion(time(2:end-idx+1), 1, 0, tau);
 
%define new measured trace
Cu = nan(1, N+1);
Cu(1)=1;
%Cu(time<=7) = 1;

for i=1:N
    if i <= idx
        dCb = photobleaching(Cu(i), alpha);
        dCu(i) = -dCb * dt;
        Cu(i+1) = Cu(i) + dCu(i);
    else
        dCb = photobleaching(Cu(i), alpha);
        dCp = diffusionDerv(time(i-idx+1), 1, 0, tau)
        dCu(i) = (-dCb + dCp) * dt;
        Cu(i+1) = Cu(i) + dCu(i);        
    end        
end

%pre-allocate vectors
Cnew = nan(1, N+1); %the corrected trace
unbFrac = nan(1, N+1); %the fractiono of unbleached fluor.
Cbl = nan(1, N+1); %the amount of bleached fluor.

Cnew(1) = Cu(1);
unbFrac(1) = 1;
Cbl(1) = 0;
   
for i=2:N

    if i<=idx
        dCB = photobleached(Cu(i-1), dt, beta, b); %this is the amount of photobleaching that occured in the bleached trace
        Cnew(i) = Cu(i) + dCB;
        Cbl(i) = Cbl(i) + dCB; %Accounting fo the change in concentration fo bleached fluorophores
        unbFrac(i)=(Cu(i))/(Cu(i)+Cbl(i)); %Calculate the new fraction of unbleached fluorophores
    elseif i>idx 
        dCB = photobleached(Cu(i-1), dt, beta, b);
        dCT = Cu(i) - Cu(i-1); %this is the total fluor. loss for the bleached trace
        dCP = dCT + dCB; %this is the amount of loss attributable to permeability
        dCP = dCP * unbFrac(i); %Correcting for the fact that a fraction of fluorophores are unbleached
        Cnew(i+1) = Cnew(i) + dCP;
        Cbl(i+1) = Cbl(i) + dCB + dCP * (1-unbFrac(i)); %Accounting fo the change in concentration fo bleached fluorophores
        unbFrac(i+1)=(Cu(i+1))/(Cu(i+1)+Cbl(i+1)); %Calculate the new fraction of unbleached fluorophores
    end
end 

%plot the 'true' and photobleached fluor. traces
figure, hold on
plot(time, C, 'Color', colorcode{5}, 'LineWidth', 1.5)
plot(time, Cu, 'Color', colorcode{7}, 'LineWidth', 1.5)
plot(time, Cnew, '--', 'Color', colorcode{7}, 'LineWidth', 1.5)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
ylim([0 Inf])


