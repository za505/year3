%diffusionAnalysis.m
%Author: Zarina Akbary
%Date: 01 May 2022
%purpose: to analyze mNeonGreen.m diffusion out of the cell

clear, close all

%color palette
okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

%% User Input
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/aggregate/';
dirsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/figures/CZI/';

basenames = {'12082021_Exp1', '04242022_Exp2', '04052022_Exp2', '04052022_Exp1'};
labels = {'Frame Rate = 2 s','Frame Rate = 3 s', 'Frame Rate = 10 s', 'Frame Rate = 20 s'};

%% Load the data and normalize
[cellTrace, bgTrace, adjTrace, times, lCell, frameRate] = loadData(dirpath, basenames);
[normTime, normTrace] = dataNormalize(times, adjTrace);

%% calculate alpha
% alpha0 = 1;
% 
% [alpha, alpha_yhat] = alphaCalc(normTime, normTrace, alpha0);
% [alphaCoef, alphaFit] = alphaLinear(frameRate, alpha);

%% calculate rho
param0 = [0, 0.01];

[rho, B, rho_yhat] = rhoCalc(normTrace,  param0);
[rhoCoef, rhoFit] = rhoLinear(frameRate, rho);

%% plot the fit of alpha
p = 0;
n1 = 1;
n2 = 4;

% figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
% for i=1:length(basenames)
%     
%     p = p + 1;
%     subplot(n1, n2, p)
%     x = normTrace{i, 1};
%     y1 = normTrace{i, 2};
%     y2 = alpha_yhat{i};
%     
%     plot(x, y1, 'Color', okabeIto{i}), hold on
%     plot(x, y2, '--k')
%     ylim([0 Inf])    
%     xlabel('Time (minutes)')
%     ylabel('Fluorescence (A.U.)')
%     title([labels{i}])
% end

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(basenames)
    
    p = p + 1;
    subplot(n1, n2, p)
    x = [0:length(normTime{i})-1];
    y1 = normTrace{i};
    y2 = rho_yhat{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, '--k')
    ylim([0 Inf])    
    xlabel('Frame Number')
    ylabel('Normalized Fluorescence (A.U.)')
    title([labels{i}])
end

%% correct the control traces
[correctedTrace]=photoCorrect(normTime, normTrace, rho1);
[correctedTrace2]=photoCorrect(normTime, normTrace2, rho1);

%% plot the control traces
p = 0;
n1 = 1;
n2 = 3;

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1 %:length(basenames)
    
    p = p + 1;
    subplot(n1, n2, p)
    x = times{i};
    y1 = cellTrace{i};
    y2 = bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(x(lysis{i}(1, 1)), '--k') %the time, not the frame
    xline(x(lysis{i}(1, 2)), '--k') %the time, not the frame
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(['Raw ' labels{i}])
    
    p = p + 1;
    subplot(n1, n2, p)
    x = normTime{i};
    y = normTrace{i};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 Inf])    
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(['Normalized ' labels{i}])
    
    p = p + 1;
    subplot(n1, n2, p)
    x = normTime{i};
    y = correctedTrace{i};
    [nr, nc] = size(y);
    ysmooth = nan(size(y));
    for j=1:nr
        ysmooth(j, :) = movingaverage(y(j, :), 3);
    end
    
    plot(x, ysmooth, 'Color', okabeIto{i})
    ylim([0 5])    
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title(['Corrected ' labels{i}])
end

%% plot the controls traces 
p = 0;
n1 = 1;
n2 = 3;

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1 %:length(basenames)
    
    p = p + 1;
    subplot(n1, n2, p)
    x = times{i};
    y1 = cellTrace{i};
    y2 = bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(x(lysis{i}(1, 1)), '--k') %the time, not the frame
    xline(x(lysis{i}(1, 2)), '--k') %the time, not the frame
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(['Raw ' labels{i}])
    
    p = p + 1;
    subplot(n1, n2, p)
    x = normTime{i};
    y = normTrace2{i};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 Inf])    
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(['Normalized ' labels{i}])
    
    p = p + 1;
    subplot(n1, n2, p)
    x = normTime{i};
    y = correctedTrace2{i};
    [nr, nc] = size(y);
    ysmooth = nan(size(y));
    for j=1:nr
        ysmooth(j, :) = movingaverage(y(j, :), 3);
    end
    
    plot(x, ysmooth, 'Color', okabeIto{i})
    ylim([0 1.5])    
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title(['Corrected ' labels{i}])
end

%% plot rho vs frame rate
dx = [];
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:4   
    x = frameRate{i}(end);
    dx = [dx x];
    y = rho{i}';
    ybar = mean(y, 2, 'omitnan');
    err = std(y, 0, 2, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')
end  
plot([0, dx], rhoFit, '--k', 'LineWidth', 1.5)

% dx = [];
% figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
% for i=1:3    
%     x = frameRate{i}(end);
%     dx = [dx x];
%     y = B{i}';
%     ybar = mean(y, 2, 'omitnan');
%     err = std(y, 0, 2, 'omitnan');
%     
%     scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
%     scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
%     errorbar(x, ybar, err, 'Color', 'black')
% end  
% plot([0, dx], rhoFit, '--k', 'LineWidth', 1.5)

%% Functions
function [cellTrace, bgTrace, adjTrace, times, lCell, frameRate] = loadData(dirpath, basenames)

    cellTrace = cell(length(basenames), 1);
    bgTrace = cell(length(basenames), 1);
    adjTrace = cell(length(basenames), 1);
    times = cell(length(basenames), 1);
    lCell = cell(length(basenames), 1);
    frameRate = cell(length(basenames), 1);

    for i=1:length(basenames)

        cd(dirpath);
        basename = basenames{i};
        datadir=dir([basename '*']);
        load(datadir(1).name, '-regexp', 'intensity$', 'time', 'lcell');

        cellTrace{i} = intensity;
        bgTrace{i} = bgintensity;
        adjTrace{i} = adjintensity;
        times{i} = time;
        lCell{i} = lcell;
        frameRate{i} = diff(time, 1, 2);
        clear intensity adjintensity bgintensity time

    end

end

function y = movingaverage(x,s)

        [mx,lx]=size(x);
        y=zeros(mx,lx);

        if s>lx
            s=round(lx/2);
        end

        for i=1:lx
                s2=floor(s/2);
            if i-s2<=0
                y(:,i)=mean(x(:,1:i+s2), 2, 'omitnan');
            elseif i+s2>lx
                y(:,i)=mean(x(:,i-s2:end), 2, 'omitnan');
            else
                y(:,i)=mean(x(:,i-s2:i+s2), 2, 'omitnan');
            end
        end
    
end

function [normTime, normTrace] = dataNormalize(times, adjTrace, imstart)
    
    %pre-allocate variables
    [nrow, ~] = size(adjTrace);
    normTime = cell(nrow, 1);
    normTrace = cell(nrow, 1);
    
    for i=1:height(adjTrace)  
        
        adjintensity = adjTrace{i};
        tme = times{i};
        
        dt = diff(tme, 1, 2);
        
        [nrow, ncol] = size(adjintensity);
        
        if nargin < 3
            imstart = find(dt == dt(end), 1, 'first');
            prelysis = imstart;
        else
            prelysis = imstart;
            lval = [prelysis + 1 : prelysis + 4];
            nval = setdiff(1:ncol, lval); %non-lysis frames
            pval = setdiff(prelysis:ncol, lval);
            fval = setdiff(prelysis:ncol-1, lval);
        
            xq = tme(lval); %the lysis times
            x = setdiff(tme, xq); %the non-lysis times
        
            for j=1:nrow
                vq = interp1(x, adjintensity(j, nval), xq);
                adjintensity(j, lval) = vq;
            end
        end
        
        
        %adjust the time vector
        normTime{i,1} = tme(prelysis:end) - tme(prelysis);
        
        %normalize the trace by the pre-lysis frame
        normTrace{i, 1} = adjintensity(:, prelysis:end) ./ adjintensity(:, prelysis);

    end
    
end

function [alpha, alpha_yhat] = alphaCalc(times, normTrace, alpha0)
    
    %pre-allocate variables
    alpha = cell(height(normTrace), 1);
    alpha_yhat = cell(height(normTrace), 1);
    
    for i=1:height(normTrace)

        normintensity = normTrace{i};
        tme = times{i};
        
        %find out if the frame rate changes
        dt = diff(tme, 1, 2);
        
        %convert negative values to NaN
        normintensity(normintensity<0) = NaN;
    
        %index traces that have at least 1 data point
        [nrow, ncol]=size(normintensity); 
        idx=find(sum(isnan(normintensity), 2) < ncol-1); %the sum of the NaNs should be 1 less than the number of time points
        normintensity=normintensity(idx, :);
        
        %remove traces without a start value
        if dt(1)~=dt(end)
            imstart = find(dt == dt(end), 1, 'first');
        else 
            imstart = 1;
        end
        
        idx = find(~isnan(normintensity(:, imstart))); 
        normintensity = normintensity(idx, :);
        
        %pre-allocate output variables
        [nrow, ~]=size(normintensity);
        alphaTemp=nan(nrow,1);    
        yhatTemp=nan(size(normintensity)); 
    
        %define exponential function
        %modelfun = @(param, t)param(1)*exp(-t./param(2));
        
        if dt(1)~=dt(end)
            idx = find(dt == dt(end), 1, 'first');
            xq = idx:ncol;
            vq = setdiff(1:ncol, xq);
            
            %fit data to the function and evaluate
            for j=1:nrow
                A = normintensity(j,idx);
                modelfun = @(alpha, t)A*exp(-t./alpha) + (1-A);
                alphaTemp(j, 1)=nlinfit(tme(xq)-tme(idx), normintensity(j,xq), modelfun, alpha0);
                yhatTemp(j, xq)=modelfun(alphaTemp(j), tme(xq)-tme(idx));
                yhatTemp(j, vq)=NaN;
            end
        else
            %fit data to the function and evaluate
            modelfun = @(alpha, t)exp(-t./alpha);
            for j=1:nrow
                alphaTemp(j, 1)=nlinfit(tme, normintensity(j,:), modelfun, alpha0);
                yhatTemp(j, :)=modelfun(alphaTemp(j), tme);
            end            
        end
   
        
        alpha{i} = alphaTemp;
        alpha_yhat{i} = yhatTemp;
        
    end    
    
end

function [alphaCoef, alphaFit] = alphaLinear(frameRate, alpha)
    
    dx = [];
    rx = [];
    xv = [];
    yv = [];

    for i=1:height(alpha)
        x = frameRate{i}(end);
        y = alpha{i}';
        dx = [dx, x];

        rx = repelem(x, length(y));
        xv = [xv, rx];
        yv = [yv, y];
    end

    alphaCoef = polyfit(xv, yv, 1);
    alphaFit = polyval(alphaCoef,[0 dx]);

end

function [rho, B, rho_yhat] = rhoCalc(normTrace, param0)
    
    %pre-allocate variables
    rho = cell(height(normTrace), 1);
    B = cell(height(normTrace), 1);
    rho_yhat = cell(height(normTrace), 1);
    
    for i=1:height(normTrace)
        
        normintensity = normTrace{i};
        
        %convert negative values to NaN
        normintensity(normintensity<0) = NaN;
    
        %index traces that have at least 1 data point
        [nrow, ncol]=size(normintensity); 
        idx=find(sum(isnan(normintensity), 2) < ncol-1); %the sum of the NaNs should be 1 less than the number of time points
        normintensity=normintensity(idx, :);
        frames = 0:ncol-1;
        
        %pre-allocate output variables
        [nrow, ~]=size(normintensity);
        rhoTemp=nan(nrow,2);    
        yhatTemp=nan(size(normintensity));
    
        %find the last time point of interest
%         di = diff(normintensity, 2, 2);
%         mdi = abs(mean(di, 1, 'omitnan'));
%         [~, fidx] = min(mdi);
        
        mint = mean(normintensity, 1, 'omitnan');
        fidx = find(mint >= 0.15, 1, 'last');
        
        %define exponential function    
        modelfun = @(param, x)(1-param(2)) * exp(-x*param(1)) + param(2);
        x = frames(1:fidx);
        
        %fit data to the function and evaluate
        for j=1:nrow
            y = normintensity(j,1:fidx);
            paramTemp(j, :)=nlinfit(x, y, modelfun, param0);
            yhatTemp(j, 1:fidx)=modelfun(paramTemp(j, :), x);
            yhatTemp(j, fidx+1:end)=NaN;
        end
                          
        rho{i} = paramTemp(:, 1);
        B{i} = paramTemp(:, 2);
        rho_yhat{i} = yhatTemp;  
        
    end  
    
end

function [rhoCoef, rhoFit] = rhoLinear(frameRate, rho)
    
    dx = [];
    rx = [];
    xv = [];
    yv = [];

    for i=1:height(rho)
        x = frameRate{i}(end);
        y = rho{i}';
        dx = [dx, x];

        rx = repelem(x, length(y));
        xv = [xv, rx];
        yv = [yv, y];
    end

    rhoCoef = polyfit(xv, yv, 1);
    rhoFit = polyval(rhoCoef,[0 dx]);

end

function [correctedTrace]=photoCorrect(normTime, normTrace, rho)
        
         correctedTrace = cell(height(normTrace), 6);
        
        for h=1:height(normTrace)
            
            normintensity = normTrace{h};
            tme = normTime{h};
            
        %pre-allocate variables
        %assume that the initial 'measured' fluorescence values and corrected
        %fluor. values will be equal. I prefer to pre-allocate with nan in case 
        %some values are missing in the raw data
        Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
        Cnew(:, 1)=normintensity(:,1);
             
        dCB=nan(height(normintensity), length(tme)-1); %change in fluor. due to photobleaching
        dCT=nan(height(normintensity), length(tme)-1); %total change in fluor.
        dCP=nan(height(normintensity), length(tme)-1); %this is the dCP, or loss attributable to permeability

        unb_frac=nan(size(normintensity)); %fraction of unbleached fluor. 
        unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point

        Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
        Cbl_exp(:, 1)=0;
        
        %this formula comes from the slope and intercept calculated for the 1.2, 2,
        %and 3 second tau vs frame rate controls 
       %dC=@(C, beta, dt, b)(C/(beta*dt+b))*dt;
       %dC=@(C, alpha, dt, b)((C*exp(-dt/(alpha*dt+b)))*(1/(alpha*dt+b))*dt);
        %dC=@(C, beta, t, dt, b)(C/(beta*t+b))*dt;
       %dC=@(C, omega, dt, b)(C*(phi*dt+b))*dt;
       
        %calculate dt (the dt between frames may vary in a single run).
        %Note: this is also the frame rate
        
        %the correction
        for n=1:height(normintensity)
           
            for i=1:length(tme)-1
                
                %dCB(n,i) = dC(normintensity(n,i), beta, tme(i), dt(i), intercept); %this is the amount of photobleaching that occured in our measured value
   
                %dCB(n,i) = dC(normintensity(n,i), phi, dt(i), intercept);
                %rho = phi + intercept;
                
                dCB(n,i) = normintensity(n,i) * rho;
                
                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value

                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                %dCP(n,i)=dCP(n,i)*unb_frac(n,i); %Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
                
            end  
            
        end
        
        correctedTrace{h, 1} = Cnew;
        correctedTrace{h, 2} = dCB;
        correctedTrace{h, 3} = dCT;
        correctedTrace{h, 4} = dCP;
        correctedTrace{h, 5} = Cbl_exp;
        correctedTrace{h, 6} = unb_frac;
        end
end

function [tau, tau_yhat] = tauCalc(times, correctedTrace, tau0)
    
    tau = cell(height(correctedTrace), 1);
    tau_yhat = cell(height(correctedTrace), 1);
    
    for i = 1:height(correctedTrace)
        x = times{i};
        y = correctedTrace{i};
          
        ybar = mean(y, 1, 'omitnan');
        A = ybar(1);
        B = ybar(end);
        
        %modelfun = @(tau, t) A * exp(-t./tau) + B;
        modelfun = @(tau, t) A * exp(log(A/B)*(1-exp(t./tau)));
        
        tau{i} = nlinfit(x, ybar, modelfun, tau0);
        tau_yhat{i} = modelfun(tau{i}, x);
    end
    
end