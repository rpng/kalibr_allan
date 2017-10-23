
% IMU information (todo: move this to the yaml file)
update_rate = 400;

% xsens mti-G-700 (page 32 and 33)
% https://www.xsens.com/download/usermanual/MTi_usermanual.pdf
gyroscope_noise_density = 0.01; %datasheet => deg/s/sqrt(Hz)
accelerometer_noise_density = 60; %datasheet => µg/sqrt(Hz)

% Convert to kalibr format
gyroscope_noise_density = pi/180*gyroscope_noise_density; %convert to => rad/s/sqrt(Hz)
accelerometer_noise_density = 9.81*accelerometer_noise_density/1000000; %convert to => m/s^2/sqrt(Hz)



%% Output results for yaml file
fprintf('\n\n=================================\n');
fprintf('rostopic: /imu0\n');
fprintf('update_rate: %.2f\n',update_rate);
fprintf('#Accelerometers\n');
fprintf('accelerometer_noise_density: %.5f.\n',accelerometer_noise_density);
fprintf('accelerometer_random_walk: %.5f.\n',0.0);
fprintf('#Gyroscopes\n');
fprintf('gyroscope_noise_density: %.6f.\n',gyroscope_noise_density);
fprintf('gyroscope_random_walk: %.6f.\n',0.0);
fprintf('=================================\n');










