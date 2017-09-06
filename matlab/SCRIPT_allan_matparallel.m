%% Initalization
close all
clear all

% Read in our toolboxes
addpath('functions')

% Our bag information
mat_path = '../data/imu.mat';

% IMU information (todo: move this to the yaml file)
update_rate = 400;

% xsens mti-G-700 (page 32 and 33)
% https://www.xsens.com/download/usermanual/MTi_usermanual.pdf
gyroscope_noise_density = 0.01; %datasheet => deg/s/sqrt(Hz)
accelerometer_noise_density = 60; %datasheet => µg/sqrt(Hz)

% Convert to kalibr format
gyroscope_noise_density = pi/180*gyroscope_noise_density; %convert to => rad/s/sqrt(Hz)
accelerometer_noise_density = 10e-6*9.81*accelerometer_noise_density; %convert to => m/s^2/sqrt(Hz)


%% Data processing
% Load the mat file (should load "data_imu" matrix)
fprintf('opening the mat file.\n')
load(mat_path);

% Load our time series information
fprintf('loading timeseries.\n')
ts_imua = timeseries(data_imu(:,2:4),data_imu(:,1));
ts_imuw = timeseries(data_imu(:,5:7),data_imu(:,1));


%% Process the timeseries data

% Find the frequency of the imu unit
% delta = mean(diff(ts_imua.Time(1:10)));
% update_rate = 1/delta;
delta = 1/update_rate;
fprintf('imu frequency of %.2f.\n',update_rate);
fprintf('sample period of %.5f.\n',delta);

% Calculate our tau range (max is half of the total measurements)
taumax = floor((length(ts_imua.Time)-1)/2);
%tau = logspace(log10(delta),log10(taumax),taumax);
tau = delta*linspace(1,taumax,500);


%% Calculate the acceleration allan deviation of the time series data!

data1.rate = update_rate;
data1.freq = ts_imua.data(:,1)';
data2.rate = update_rate;
data2.freq = ts_imua.data(:,2)';
data3.rate = update_rate;
data3.freq = ts_imua.data(:,3)';

fprintf('calculating acceleration allan deviation.\n');
tic;
j1 = batch(@allan,1,{data1,tau});
j2 = batch(@allan,1,{data2,tau});
j3 = batch(@allan,1,{data3,tau});

% Wait for the jobs to finish
wait(j1)
wait(j2)
wait(j3)

% Get results into a cell array
r1 = fetchOutputs(j1);
r2 = fetchOutputs(j2);
r3 = fetchOutputs(j3);
results_ax = r1{1};
results_ay = r2{1};
results_az = r3{1};
toc;

% Plot the results on a figure
figure(1);
loglog(results_ax.tau1, sqrt(results_ax.sig2)); hold on;
loglog(results_ay.tau1, sqrt(results_ay.sig2)); hold on;
loglog(results_az.tau1, sqrt(results_az.sig2)); hold on;
grid on;
xlabel('\tau [sec]');
ylabel('Normal Allan Deviation [m/s^2]');
legend('x-acceleration','y-acceleration','z-acceleration');

pause(5)


%% Calculate the angular allan deviation of the time series data!

data1.rate = update_rate;
data1.freq = ts_imuw.data(:,1)';
data2.rate = update_rate;
data2.freq = ts_imuw.data(:,2)';
data3.rate = update_rate;
data3.freq = ts_imuw.data(:,3)';

fprintf('calculating angular allan deviation.\n');
tic;
j1 = batch(@allan,1,{data1,tau});
j2 = batch(@allan,1,{data2,tau});
j3 = batch(@allan,1,{data3,tau});

% Wait for the jobs to finish
wait(j1)
wait(j2)
wait(j3)

% Get results into a cell array
r1 = fetchOutputs(j1);
r2 = fetchOutputs(j2);
r3 = fetchOutputs(j3);
results_wx = r1{1};
results_wy = r2{1};
results_wz = r3{1};
toc;

% Plot the results on a figure
figure(2);
loglog(results_wx.tau1, sqrt(results_wx.sig2)); hold on;
loglog(results_wy.tau1, sqrt(results_wy.sig2)); hold on;
loglog(results_wz.tau1, sqrt(results_wz.sig2)); hold on;
grid on;
xlabel('\tau [sec]');
ylabel('Normal Allan Deviation [rad/s]');
legend('x-angular','y-angular','z-angular');



%% Output results for yaml file
fprintf('\n=================================\n');

fprintf('rostopic: /imu0\n');
fprintf('update_rate: %.2f\n',freq);
fprintf('\n#Accelerometers\n');
fprintf('accelerometer_noise_density: %.5f.\n',accelerometer_noise_density);
fprintf('accelerometer_random_walk: %.5f.\n',0.0);
fprintf('\n#Gyroscopes\n');
fprintf('gyroscope_noise_density: %.5f.\n',gyroscope_noise_density);
fprintf('gyroscope_random_walk: %.5f.\n',0.0);






