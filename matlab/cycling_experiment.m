s = serial('COM3');
fopen(s);

line = fgetl(s);

duty_cycles = [0:0.1:1 0.9:-0.1:0];
m = length(values);
for i = 1:m
    % Set new duty cycle
    d = min(max(duty_cycles(i),0),1);    
    fwrite(s,uint8(d * 255));
    for j = 1:n
        line = fgetl(s);
        % Save data
    end
    % Take image    
end