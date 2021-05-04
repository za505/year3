%%PlotScrapBook.m
%Purpose: to save code for one-off plots I've made

%let's plot the average ratio too
figure, hold on
title('Intensity/Background vs Time')
xlabel('Time (s)')
ylabel('Intensity/Background')
ylim([0 1])
fig2pretty
%yline(1, '--k')
for x=1:length(xswitch)
   xline(xswitch(x), '--k') 
end
ciplot(avgRatio - stdRatio, avgRatio + stdRatio, time, [0.75 0.75 1])
plot(time, avgRatio, '-r') 