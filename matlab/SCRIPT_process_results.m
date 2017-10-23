%% Initalization
close all
clear all

% Read in our toolboxes
addpath('functions/allan_v3')

% Our bag information
mat_path = '../data/results.mat';

% IMU information (todo: move this to the yaml file)
update_rate = 400;


%% Data processing
% Load the mat file (should load "data_imu" matrix)
fprintf('opening the mat file.\n')
load(mat_path);

% Remove NANs
% results_ax.tau1(~any(~isnan(results_ax.sig2), 2),:)=[];
% results_ax.sig2(~any(~isnan(results_ax.sig2), 2),:)=[];

% Calculate norms
results_aavg.tau1 = results_ax.tau1;
results_aavg.sig2 = mean([results_ax.sig2;results_ay.sig2;results_az.sig2]);


%% Plot the results on a figure
fh1 = figure(1);
loglog(results_ax.tau1, sqrt(results_ax.sig2)); hold on;
loglog(results_ay.tau1, sqrt(results_ay.sig2)); hold on;
loglog(results_az.tau1, sqrt(results_az.sig2)); hold on;
loglog(results_aavg.tau1, sqrt(results_aavg.sig2)); hold on;
grid on;
%axis equal;
xlabel('\tau [sec]');
ylabel('Normal Allan Deviation [m/s^2]');
legend('x-acceleration','y-acceleration','z-acceleration','average');


% =======================================================================
% Find location of tau=1 by finding where difference is near zero
tauref = 1;
taudiff = abs(results_ax.tau1-tauref);
tauid = find(taudiff == min(taudiff));
%fprintf('taua = %.2f | taua-id = %d\n',results_ax.tau1(tauid),tauid);


% =======================================================================
% We will fit our "white-noise" line where tau is 1
windowsize = 20;
window = tauid-windowsize:tauid+windowsize;
x = results_ax.tau1(window);
y = sqrt(results_aavg.sig2(window));
nanx = isfinite(y);
% Fit to log-log scale
coeffs = polyfit(log(x(nanx)), log(y(nanx)), 1);
%fprintf('slope = %.4f | y-intercept = %.4f\n',coeffs(1),coeffs(2));
% Convert from logarithmic scale to linear scale and plot
h_fit = results_ax.tau1.^coeffs(1).*exp(coeffs(2));
plot(results_ax.tau1, h_fit,'Color','r','LineWidth',2); hold on;


% =======================================================================
% Get our "white-noise" value at the fitted line, and from the raw data
fprintf('accelerometer_noise_density = %.6f [raw data]\n',sqrt(results_aavg.sig2(tauid)));
fprintf('accelerometer_noise_density = %.6f [fitted line]\n',h_fit(tauid));



% =======================================================================
% Next we should fit a -1/2 line to the average







%% Plot the results on a figure
fh2 = figure(2);
loglog(results_wx.tau1, sqrt(results_wx.sig2)); hold on;
loglog(results_wy.tau1, sqrt(results_wy.sig2)); hold on;
loglog(results_wz.tau1, sqrt(results_wz.sig2)); hold on;
grid on;
xlabel('\tau [sec]');
ylabel('Normal Allan Deviation [rad/s]');
legend('x-angular','y-angular','z-angular');


% Find location of tau=1 by finding where difference is near zero
tauref = 1;
taudiff = abs(results_wx.tau1-tauref);
tauid = find(taudiff == min(taudiff));
%fprintf('taug = %.2f | taug-id = %d\n',results_wx.tau1(tauid),tauid);
siggs = [sqrt(results_wx.sig2(tauid)),sqrt(results_wy.sig2(tauid)),sqrt(results_wz.sig2(tauid))];
fprintf('= sigg_x: %.6f\n',siggs(1));
fprintf('= sigg_y: %.6f\n',siggs(2));
fprintf('= sigg_z: %.6f\n',siggs(3));
fprintf('gyroscope_noise_density = %.6f\n\n',mean(siggs));



% Save to file
%print(fh1,'-dpng','-r500','../data/allan_linear_acceleration.png')
%print(fh2,'-dpng','-r500','../data/allan_angular_velocity.png')
