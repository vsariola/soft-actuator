function compactdaq_test

    s = daq.createSession('ni');
    s.Rate = 50000;
    s.DurationInSeconds = 3;
    addChannel('ai0');
    addChannel('ai1');
    addChannel('ai2');  
    
    lh = addlistener(s,'DataAvailable',@plotData); 
 
    function plotData(src,event)
        plot(event.TimeStamps,event.Data)
    end
    startBackground(s);
    pause;

    function ch = addChannel(channel_name)
        ch = addAnalogInputChannel(s,'cDAQ1Mod1',channel_name,'Bridge');
        ch.BridgeMode = 'Half';
        ch.NominalBridgeResistance = 120;        
    end
end

