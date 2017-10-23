%% Initalization
close all
clear all

% Read in our toolboxes
addpath('functions/allan_v3')

% Our bag information
bag_path = '../data/imu.bag';
imu_topic = '/sensor/xsens/sensor/imu_raw';

% Bag start and end time
bagstart = 0;
bagend = 30;


%% Data processing
% Open ros bag, and select topics we want
fprintf('opening the ros bag.\n')
filepath = fullfile(bag_path);
bag = rosbag(filepath);

% Select the topics we should insert
fprintf('selecting topics.\n')
bagselect = select(bag, 'Time', [bag.StartTime+bagstart bag.StartTime+bagend], 'Topic', imu_topic);

% Load our time series information
fprintf('loading timeseries.\n')
ts_imuw = timeseries(bagselect, 'AngularVelocity.X', 'AngularVelocity.Y', 'AngularVelocity.Z');
ts_imua = timeseries(bagselect, 'LinearAcceleration.X', 'LinearAcceleration.Y', 'LinearAcceleration.Z');

% Change the start time so that everything is relative to the first one
%ts_imuw.TimeInfo.StartDate = datestr(datetime(ts_imuw.Time(1), 'ConvertFrom', 'posixtime'));
%ts_imua.TimeInfo.StartDate = datestr(datetime(ts_imua.Time(1), 'ConvertFrom', 'posixtime'));


%% Process the timeseries data

% Find the frequency of the imu unit
% delta = mean(diff(ts_imua.Time(1:10)));
% freq = 1/delta;
freq = 400;
delta = 1/freq;
fprintf('imu frequency of %.2f.\n',freq);
fprintf('sample period of %.5f.\n',delta);

% Set what our frequency is for the system
data.rate = freq;

% Calculate our tau range (max is half of the total measurements)
taumax = (length(ts_imua.Time)-1)/2;
tau = delta*linspace(1,taumax,5000);


% Calculate the allan deviation of the time series data!
fprintf('calculating acceleration x allan deviation.\n');
data.freq = ts_imua.data(:,1)';
[results_ax] = allan(data, tau);

fprintf('calculating acceleration y allan deviation.\n');
data.freq = ts_imua.data(:,2)';
[results_ay] = allan(data, tau);

fprintf('calculating acceleration y allan deviation.\n');
data.freq = ts_imua.data(:,3)';
[results_az] = allan(data, tau);


%% Plot the results on a figure
figure;
loglog(results_ax.tau1, sqrt(results_ax.sig2)); hold on;
loglog(results_ay.tau1, sqrt(results_ay.sig2)); hold on;
loglog(results_az.tau1, sqrt(results_az.sig2)); hold on;

grid on;
xlabel('\tau [sec]');
ylabel('Normal Allan Deviation [m/s^2]');
legend('x-accel','y-accel','z-accel');









