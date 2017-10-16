curvature = analyzeVideo('D:\13.10.2017 tests for microactuators lab\MVI_4352.mov');
load('D:\13.10.2017 tests for microactuators lab\cycling_data.mat');
pressure = data(:,2);
[pa,ca] = alignData(pressure,curvature);
plot(pa,ca);