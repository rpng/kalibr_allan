%% Initalization
close all
clear all

% Read in our toolboxes
addpath('functions/allan_v3')

% Our bag information
%mat_path = '../data/imu_mtig700.mat';
%mat_path = '../data/imu_tango.mat';
mat_path = '../data/imu_visensor.mat';

% IMU information (todo: move this to the yaml file)
%update_rate = 400;
%update_rate = 100;
update_rate = 200;


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
tau = delta*logspace(log10(delta),log10(taumax),2000);
%tau = delta*linspace(1,taumax,1000);


%% Calculate the acceleration allan deviation of the time series data!

data1.rate = update_rate;
data1.freq = ts_imua.data(:,1)';
data2.rate = update_rate;
data2.freq = ts_imua.data(:,2)';
data3.rate = update_rate;
data3.freq = ts_imua.data(:,3)';
data4.rate = update_rate;
data4.freq = ts_imuw.data(:,1)';
data5.rate = update_rate;
data5.freq = ts_imuw.data(:,2)';
data6.rate = update_rate;
data6.freq = ts_imuw.data(:,3)';

fprintf('calculating allan deviation.\n');
tic;
cluster = parcluster();
j1 = batch(cluster,@allan,1,{data1,tau});
j2 = batch(cluster,@allan,1,{data2,tau});
j3 = batch(cluster,@allan,1,{data3,tau});
j4 = batch(cluster,@allan,1,{data4,tau});
j5 = batch(cluster,@allan,1,{data5,tau});
j6 = batch(cluster,@allan,1,{data6,tau});

% Wait for the jobs to finish
wait(j1)
wait(j2)
wait(j3)
wait(j4)
wait(j5)
wait(j6)

% Get results into a cell array
r1 = fetchOutputs(j1);
r2 = fetchOutputs(j2);
r3 = fetchOutputs(j3);
r4 = fetchOutputs(j4);
r5 = fetchOutputs(j5);
r6 = fetchOutputs(j6);
results_ax = r1{1};
results_ay = r2{1};
results_az = r3{1};
results_wx = r4{1};
results_wy = r5{1};
results_wz = r6{1};

% Finally cleanup
delete(j1)
delete(j2)
delete(j3)
delete(j4)
delete(j5)
delete(j6)
toc

%% Save workspace
filename = ['results_',datestr(now,30),'.mat'];
fprintf('saving to: %s\n',filename);
save(['../data/',filename],'update_rate','ts_imua','ts_imuw','tau','taumax','results_ax','results_ay','results_az','results_wx','results_wy','results_wz')
fprintf('done saving!\n');




