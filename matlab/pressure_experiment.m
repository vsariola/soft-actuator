n = 1000;

delete(instrfindall); % delete all old objects from memory

data = zeros(n,6);

s = serial('COM10');
fopen(s);

fgetl(s);
tic;
for i = 1:n
    line = fgetl(s);
    linevals = cell2mat(textscan(line,'%f'));
    data(i,:) = [linevals' toc];
end
fclose(s);

save('pressure_data.mat','data');