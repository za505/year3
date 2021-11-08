%% Modeling Photobleach Correction
%Rico's de-bugging
clear, close all
%% Modeling Photobleach Correction, 1 minute frame rate

%set variables
N = 100; %number of minutes
dt1 = 1; %increment in minutes, frame rate = 1 minute
T1 = [0:dt1:N]; %time vector in minutes


%I hypothesize that fluorescence decreases according to the following
%equation: dC/dt = -beta * C + gamma * C, where C is the measured
%fluorescence, beta is the rate of photobleaching, and gamma is the
%permeability coefficient
gamma = 0.01; %permeability coefficient (I made gamma positive)
beta1 = 1/16.4592; %photobleaching coefficient, beta=1/tau
beta2 = 1/8.7189; %frame rate = 30 s
beta3 = 1/4.3229; %frame rate = 15 s

dt2=0.5; %30 s
dt3=0.25; %15 s

%pre-allocte 'measured' fluorescence variables
%assume that our initial measured fluorescence is equal to the initial 'no diffusion' fluorescence

Cbl=zeros(size(T1));%Normalized concentration of bleached fluorophores
Cunb=zeros(size(T1));%Normalized concentration of unbleached fluorophores
Ctot=zeros(size(T1));%Normalized concentration of total fluorophores
Cunb(1)=1;%Initial concentration of unbleached fluorophores is 1
Cbl(1)=0;%Initial concentration of unbleached fluorophores is 0
Ctot(1)=1;

%using the following Newton Method, I expect to generate an exponential
%decay curve
for i=1:length(T1)-1
    dCunb= -beta1 * Cunb(i) - gamma * Cunb(i);
    dCbl=beta1 * Cunb(i) - gamma * Cbl(i);
    
    dCunb=dCunb*dt1;
    dCbl=dCbl*dt1;
    
    Cunb(i+1)=Cunb(i)+dCunb;
    Cbl(i+1)=Cbl(i)+dCbl;

    Ctot(i+1)=Cunb(i+1)+Cbl(i+1);
   
end

figure, hold on
plot(T1, Cunb, '-r') %this should be a decay
plot(T1, Cbl, '-b') %this should increase exponentially? 
plot(T1, Ctot, '-g') %this should also be an exponential decay, but less steep than Cunb
legend({'Unbleached', 'Bleached', 'Total'})
%at some point, all the fluor. are bleached and no unbleached are left, but
%that does not mean there are no fluor. in the cell

%what about what the traces look like with only gamma? (aka corrected for
%photobleaching)
M1 = zeros(size(T1));
M1(1)=1;

%using the following Newton Method, we expect to generate exponential
%decay. This is what we expect to see after correcting for photobleaching,
%when there is only the effect of the permeability coeffecient 
for i=1:length(T1)-1
    m = - gamma * M1(i); 
    dM1 = m*dt1; 
    M1(i+1) = M1(i) + dM1;
end

figure, hold on
plot(T1, M1, '-r')
plot(T1, Ctot, '--g')
legend({'M1', 'Ctot'}) %the traces completely overlap

%now that I'm able to generate exponentials, can I work in the reverse
%direction and get M and F using alpha?

%first, plot tau vs time and fit to a linear. The slope is alpha (dtau/dt,
%so 1/alpha = dbeta/dt?) I am not sure why this slope is useful or what the
%units are or exactly how it stands in place for beta?
% linearCoef = polyfit([dt3, dt2, dt1],[1/beta3, 1/beta2, 1/beta1],1);
% alpha = linearCoef(1); 

alpha=1/beta1/dt1;%Here I just used the alpha that would result from just using beta1 for my measured time constant
%this would be the same as tau though??

%assume that the initial 'measured' fluorescence values and corrected
%fluor. values will be equal
Cnew1=nan(1, length(T1));%Corrected concentration of fluorophores
Cnew1(1)=1;

%this is the dCP, or loss attributable to permeability
dCP1=nan(1, length(T1)-1);

unb_frac=zeros(size(Cunb));
unb_frac(1)=1;%Fraction of fluorophores that are unbleached

Cbl_exp=zeros(size(T1));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores

%this correction should work regardless of the frame rate
for i=1:length(T1)-1
   
    dCB = Cunb(i)/(alpha); %this is the amount of photobleaching that occured in our measured value
    dCT = Cunb(i+1) - Cunb(i); %this is the total fluor. loss for the measured value
    dCP1(i) = dCT + dCB; %this is the amount of loss attributable to permeability
    
    
    dCP1(i)=dCP1(i)/unb_frac(i);%Correcting for the fact that a fraction of fluorophores are unbleached
    
    Cnew1(i+1)=Cnew1(i)+dCP1(i);
    
    Cbl_exp(i+1)=Cbl_exp(i)+dCB+dCP1(i)*(1-unb_frac(i));%Accounting fo the change in concentration fo bleached fluorophores
     
    unb_frac(i+1)=(Cunb(i+1))/(Cunb(i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores
    
end

%plot these values. The Cnew1 vector and M vector should look the same, but
%they don't. What went wrong with the correction? 
figure, hold on
plot(T1, Cnew1, '-b')
plot(T1, M1, '--m')
plot(T1,Cunb,'.k')
legend({'Corrected Conc.','Actual Conc.','Measured Conc.'})