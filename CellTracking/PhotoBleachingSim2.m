%% Modeling Photobleach Correction
%Authors: Enrique Rojas and Zarina Akbary
%Date: 11/07/2021
%Purpose: to model the photobleach correction

clear, close all

% %Let's add some noise to our measured fluor. to see how robust the
% %adjustment is
% Cunb = Cunb + randn(1, length(Cunb));
%% Inputs

colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
%% Modeling Photobleach Correction
%bounds: 0 <= gamma <= 1
%        0 < dt

%inputs: number of time points, dt, and gamma (permeability coefficient)
%outputs: 
% T = time vector
% Cbl = concentration of bleached fluor.
% Cunb = concentration of unbleached fluor. or measured fluor.
% Ctot = total concentration of intracellular fluor.
% M = modeled fluor concentration using the Newton Method, exp(-t*gamma)
% beta = photobleaching coeff. calculated as a func of dt
% F = modeled fluor concentration using a MATLAB function, exp(-t*gamma)

%inputs: N (number of time steps), dt, gamma
%LB control, frame rate = 1 min
[T1, Cbl1, Cunb1, Ctot1, M1, beta1, F1] = simRun(100, 1, 1);

%LB control, frame rate = 30 s
[T2, Cbl2, Cunb2, Ctot2, M2, beta2, F2] = simRun(100, 0.5, 0.1);

%LB control, frame rate = 15 s
[T3, Cbl3, Cunb3, Ctot3, M3, beta3, F3] = simRun(100, 0.25, 0.05);

%PBS incubation, frame rate = 1 min
[T4, Cbl4, Cunb4, Ctot4, M4, beta4, F4] = simRun(100, 1, 0.03);

%% how much does the 'actual' differ from C total? how much do the modeled values differ from each other?
if check1==1
    figure, hold on
    plot(T1, M1, 'Color', colorcode{1}, 'LineStyle', '-') 
    plot(T1, Ctot1, 'Color', colorcode{1}, 'LineStyle', '--') 
    plot(T2, M2, 'Color', colorcode{2}, 'LineStyle', '-') 
    plot(T2, Ctot2, 'Color', colorcode{2}, 'LineStyle', '--') 
    plot(T3, M3, 'Color', colorcode{3}, 'LineStyle', '-') 
    plot(T3, Ctot3, 'Color', colorcode{3}, 'LineStyle', '--') 
    plot(T4, M4, 'Color', colorcode{4}, 'LineStyle', '-') 
    plot(T4, Ctot4, 'Color', colorcode{4}, 'LineStyle', '--') 
    title('how much does the actual (M) differ from C total?')
% legend({'LB, frame rate = 1 min, actual', 'LB, frame rate = 1 min, measured', 'LB, frame rate = 30 s, actual', 'LB, frame rate = 30 s, measured', 'LB, frame rate = 15 s, actual', 'LB, frame rate = 15 s, measured', 'PBS, frame rate = 1 min, actual', 'PBS, frame rate = 1 min, measured'})
end
% Observation: they do not differ, which is what we expect

%% how much does the 'actual' differ from the measured?
if check2==1
    figure, hold on
    plot(T1, M1, 'Color', colorcode{1}, 'LineStyle', '-') 
    plot(T1, Cunb1, 'Color', colorcode{1}, 'LineStyle', '--') 
    plot(T2, M2, 'Color', colorcode{2}, 'LineStyle', '-') 
    plot(T2, Cunb2, 'Color', colorcode{2}, 'LineStyle', '--') 
    plot(T3, M3, 'Color', colorcode{3}, 'LineStyle', '-') 
    plot(T3, Cunb3, 'Color', colorcode{3}, 'LineStyle', '--') 
    plot(T4, M4, 'Color', colorcode{4}, 'LineStyle', '-') 
    plot(T4, Cunb4, 'Color', colorcode{4}, 'LineStyle', '--') 
end

%Observation: For gamma = 0.8, the fluorphores leak out of the cell faster
%than they can be bleached. In this case, gamma is greater than beta, which
%may indicate why correction would be difficult, or perhaps even
%unnecessary

%% Can we obtain the 'actual' fluorescence value given only the measured values?
[Cnew1, dCB1, dCU1, dCP1, tau1, gamma1] = simReverse(Cunb1, T1, 1);
[Cnew2, dCB2, dCU2, dCP2, tau2, gamma2] = simReverse(Cunb2, T2, 0.5);
[Cnew3, dCB3, dCU3, dCP3, tau3, gamma3] = simReverse(Cunb3, T3, 0.25);
[Cnew4, dCB4, dCU4, dCP4, tau4, gamma4] = simReverse(Cunb4, T4, 1);

%% how close are the photocorrected values and the 'actual' fluor?
if check3==1
    figure, hold on
    plot(T1, M1, 'Color', colorcode{1}, 'LineStyle', 'none', 'Marker', 'o') 
    plot(T1, Cnew1, 'Color', colorcode{1}, 'LineStyle', '-') 
    plot(T2, M2, 'Color', colorcode{2}, 'LineStyle', 'none', 'Marker', 'o') 
    plot(T2, Cnew2, 'Color', colorcode{2}, 'LineStyle', '-') 
    plot(T3, M3, 'Color', colorcode{3}, 'LineStyle', 'none', 'Marker', 'o') 
    plot(T3, Cnew3, 'Color', colorcode{3}, 'LineStyle', '-') 
    plot(T4, M4, 'Color', colorcode{4}, 'LineStyle', 'none', 'Marker', 'o') 
    plot(T4, Cnew4, 'Color', colorcode{4}, 'LineStyle', '-') 
end

%% Does F = exp(time*-gamma) equal M?
if check4==1
    figure, hold on
    plot(T1, M1, 'Color', colorcode{1}, 'LineStyle', 'none', 'Marker', 'o') 
    plot(T1, F1, 'Color', colorcode{1}, 'LineStyle', '-') 
    plot(T2, M2, 'Color', colorcode{2}, 'LineStyle', 'none', 'Marker', 'o') 
    plot(T2, F2, 'Color', colorcode{2}, 'LineStyle', '-') 
    plot(T3, M3, 'Color', colorcode{3}, 'LineStyle', 'none', 'Marker', 'o') 
    plot(T3, F3, 'Color', colorcode{3}, 'LineStyle', '-') 
    plot(T4, M4, 'Color', colorcode{4}, 'LineStyle', 'none', 'Marker', 'o') 
    plot(T4, F4, 'Color', colorcode{4}, 'LineStyle', '-') 
    title('how much do the modeled values differ from each other?')
end

%% Functions
function [T, Cbl, Cunb, Ctot, M, beta, F] = simRun(N, dt, gamma)

% N = number of minutes
% dt = increment in minutes, frame rate 
% gamma = permeability coefficient

T = [0:dt:N]; %time vector in minutes
%tau = (16.0806)*dt + 0.4234;
tau = (30)*dt;
beta = 1/tau;

%I hypothesize that measured fluorescence decreases according to the following
%equation: dC/dt = -beta * C + gamma * C, where C is the measured
%fluorescence, beta is the rate of photobleaching, and gamma is the
%permeability coefficient. However, the true decrease in fluorophore
%concentration only varies by gamma. 

%pre-allocte 'measured' fluorescence variables
%assume that our initial measured fluorescence is equal to the initial 'no diffusion' fluorescence

Cbl=nan(size(T));%Normalized concentration of bleached fluorophores
Cunb=nan(size(T));%Normalized concentration of unbleached fluorophores
Ctot=nan(size(T));%Normalized concentration of total fluorophores
Cunb(1)=1;%Initial concentration of unbleached fluorophores is 1
Cbl(1)=0;%Initial concentration of unbleached fluorophores is 0
Ctot(1)=1;

%using the following Newton Method, I expect to generate an exponential
%decay curve
for i=1:length(T)-1
    dCunb= -beta * Cunb(i) - gamma * Cunb(i);
    dCbl= beta * Cunb(i) - gamma * Cbl(i);
    
    dCunb=dCunb*dt;
    dCbl=dCbl*dt;
    
    Cunb(i+1)=Cunb(i)+dCunb;
    Cbl(i+1)=Cbl(i)+dCbl;

    Ctot(i+1)=Cunb(i+1)+Cbl(i+1);
   
end

%what about what the traces look like with only gamma? (aka corrected for
%photobleaching)
M = zeros(size(T));
M(1)=1;

%using the following Newton Method, we expect to generate exponential
%decay. This is what we expect to see after correcting for photobleaching,
%when there is only the effect of the permeability coeffecient 
for i=1:length(T)-1
    m = - gamma * M(i); 
    dM = m*dt; 
    M(i+1) = M(i) + dM;
end

%can I get the same fluor. values as M using the equation C =
%exp(time.*-gamma)?
modelExp=@(gamma,time)exp(time.*-gamma);
F = modelExp(gamma, T);

end

function [Cnew, dCB, dCU, dCP, tau, gamma] = simReverse(Cunb, T, dt)

    %first, let's fit the measured fluor. values to an exponential to calculate
    %tau

    modelfun1=@(tau,x)exp(-x./tau);
    tau0=1;

    tau=nlinfit(T, Cunb, modelfun1, tau0);


%can I work in the reverse direction and get M using alpha (the slope of tau vs frame rate)?
%alpha=1/beta/dt; %the proportion of photobleaching
beta = 1/((16.0806)*dt + 0.4234);
alpha=1/beta/dt;

%assume that the initial 'measured' fluorescence values and corrected
%fluor. values will be equal
Cnew=nan(size(T));%Corrected concentration of fluorophores
Cnew(:, 1)=Cunb(:, 1);

%this is the dCP, or loss attributable to permeability
dCB=nan(height(Cunb), length(T)-1);
dCU=nan(height(Cunb), length(T)-1);
dCP=nan(height(Cunb), length(T)-1);

unb_frac=zeros(size(T));
unb_frac(:,1)=1;%Fraction of fluorophores that are unbleached

Cbl_exp=nan(size(T));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
Cbl_exp(:,1)=0;

if tau>=1
%this correction should work regardless of the frame rate
for i=1:length(T)-1
   
    dCB(i) = Cunb(i)/alpha; %this is the amount of photobleaching that occured in our measured value
    dCU(i) = Cunb(i+1) - Cunb(i); %this is the total fluor. loss for the measured value
    dCP(i) = dCU(i) + dCB(i); %this is the amount of loss attributable to permeability
    
    
    dCP(i)=dCP(i)/unb_frac(i);%Correcting for the fact that a fraction of fluorophores are unbleached
    
    Cnew(i+1)=Cnew(i)+dCP(i);
    
    Cbl_exp(i+1)=Cbl_exp(i)+dCB(i)+dCP(i)*(1-unb_frac(i));%Accounting fo the change in concentration fo bleached fluorophores
     
    unb_frac(i+1)=(Cunb(i+1))/(Cunb(i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores
    
end

elseif tau<1
    
    %this correction should work regardless of the frame rate
for i=1:length(T)-1
   
    dCB(i) = ((Cunb(i)+Cunb(i+1))/2)/alpha; %this is the amount of photobleaching that occured in our measured value
    dCU(i) = Cunb(i+1) - Cunb(i); %this is the total fluor. loss for the measured value
    dCP(i) = dCU(i) + dCB(i); %this is the amount of loss attributable to permeability
    
    
    dCP(i)=dCP(i)/unb_frac(i);%Correcting for the fact that a fraction of fluorophores are unbleached
    
    Cnew(i+1)=Cnew(i)+dCP(i);
    
    Cbl_exp(i+1)=Cbl_exp(i)+dCB(i)+dCP(i)*(1-unb_frac(i));%Accounting fo the change in concentration fo bleached fluorophores
     
    unb_frac(i+1)=(Cunb(i+1))/(Cunb(i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores
    
end

end

    modelfun=@(gamma,x)exp(-x.*gamma);
    gamma0=1;
    gamma=nlinfit(T, Cnew, modelfun, gamma0);

end