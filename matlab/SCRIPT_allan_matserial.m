%% Initalization
close all
clear all

% Read in our toolboxes
addpath('functions/allan_v3')

% Our bag information
mat_path = '../data/imu.mat';


% IMU information (todo: move this to the yaml file)
update_rate = 400;


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
fprintf('imu frequency of %.2f.\n',freq);
fprintf('sample period of %.5f.\n',delta);

% Set what our frequency is for the system
data.rate = update_rate;

% Calculate our tau range (max is half of the total measurements)
taumax = floor((length(ts_imua.Time)-1)/2);
%tau = logspace(log10(delta),log10(taumax),taumax);
tau = linspace(1,taumax,taumax);


%% Calculate the acceleration allan deviation of the time series data!
fprintf('calculating acceleration x allan deviation.\n');
data.freq = ts_imua.data(:,1)';
[results_ax] = allan(data, tau);

fprintf('calculating acceleration y allan deviation.\n');
data.freq = ts_imua.data(:,2)';
[results_ay] = allan(data, tau);

fprintf('calculating acceleration y allan deviation.\n');
data.freq = ts_imua.data(:,3)';
[results_az] = allan(data, tau);

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
fprintf('calculating angular x allan deviation.\n');
data.freq = ts_imuw.data(:,1)';
[results_wx] = allan(data, tau);

fprintf('calculating angular y allan deviation.\n');
data.freq = ts_imuw.data(:,2)';
[results_wy] = allan(data, tau);

fprintf('calculating angular y allan deviation.\n');
data.freq = ts_imuw.data(:,3)';
[results_wz] = allan(data, tau);



% Plot the results on a figure
figure(2);
loglog(results_wx.tau1, sqrt(results_wx.sig2)); hold on;
loglog(results_wy.tau1, sqrt(results_wy.sig2)); hold on;
loglog(results_wz.tau1, sqrt(results_wz.sig2)); hold on;
grid on;
xlabel('\tau [sec]');
ylabel('Normal Allan Deviation [rad/s]');
legend('x-angular','y-angular','z-angular');








