function [fh1, sigmaw, sigmab] = gen_chart(tau,sigx,sigy,sigz,titlestr,name,unit0,unit1,unit2)

% Calculate our average sigma values
sigavg = mean([sigx;sigy;sigz]);

% =======================================================================
% Plot the results on a figure
fh1 = figure;
loglog(tau, sigx); hold on;
loglog(tau, sigy); hold on;
loglog(tau, sigz); hold on;
loglog(tau, sigavg); hold on;
grid on;
title([titlestr,' ',name]);
xlabel('\tau [sec]');
ylabel(['Normal Allan Deviation [',unit0,']']);
legend(['x-',name],['y-',name],['z-',name],'average','Location','southeast');

% =======================================================================
% Find location of tau=1 by finding where difference is near zero
tauref = 1;
taudiff = abs(tau-tauref);
tauid1 = find(taudiff == min(taudiff),1);
fprintf('tau = %.2f | tauid1 = %d\n',tau(tauid1),tauid1);

% We will fit our "white-noise" line where tau is 1
windowsize = 50;
%window = tauid1-windowsize:tauid1+windowsize;
minid = find(sigavg == min(sigavg)); % find where the min is
window = 1:minid-windowsize; % go from start to min

% Get our x and y values
x = tau(window);
y = sigavg(window);
nanx = isfinite(y);

% Fit to log-log scale (slope of -1/2)
%coeffs = polyfit(log(x(nanx)), log(y(nanx)), 1);
coeffs(1)= -0.5;
intcs = log(y(nanx)./x(nanx).^coeffs(1));
coeffs(2) = mean(intcs);
fprintf('h_fit1 slope = %.4f | y-intercept = %.4f\n',coeffs(1),coeffs(2));

% Convert from logarithmic scale to linear scale and plot
h_fit1 = tau.^coeffs(1).*exp(coeffs(2));
pltlin = plot(tau, h_fit1,'Color','r','LineWidth',2); hold on;
pltlin.Color(4) = 0.5; % make it semi-transparent


% =======================================================================
% Next we should fit a 1/2 line to the bias side of the allan plot
minid = find(sigavg == min(sigavg)); % find where the min is
window = minid+windowsize:length(tau); % go from min to the end

% Get our x and y values
x = tau(window);
y = sigavg(window);
nanx = isfinite(y);

% Calculate the intercept given the slope of +1/2
coeffs(1)= 0.5;
intcs = log(y(nanx)./x(nanx).^coeffs(1));
coeffs(2) = mean(intcs);
fprintf('h_fit2 slope = %.4f | y-intercept = %.4f\n',0.5,mean(intcs));

% Convert from logarithmic scale to linear scale and plot
h_fit2 = tau.^coeffs(1).*exp(coeffs(2));
pltlin = plot(tau, h_fit2,'Color','b','LineWidth',2); hold on;
pltlin.Color(4) = 0.5; % make it semi-transparent


% Get what our bias should be (should be at a tau of 3)
tauref = 3;
taudiff = abs(tau-tauref);
tauid2 = find(taudiff == min(taudiff),1);
fprintf('tau = %.2f | tauid2 = %d\n',tau(tauid2),tauid2);


% =======================================================================
% Record the values for the output
sigmaw = h_fit1(tauid1);
sigmab = h_fit2(tauid2);

% Plot the values
str1 = sprintf('\\sigma = %.6f %s',h_fit1(tauid1),unit1);
str2 = sprintf('\\sigma_{b} = %.6f %s',h_fit2(tauid2),unit2);
text(.1,.15,{str1,str2},'Units','normalized'); %'FontWeight','bold'


end

