function x = toNumericIfPossible(x)
    if (ischar(x))
        foo = str2double(x);
        if (~isempty(foo) && ~isnan(foo))
            x = foo;
        end
    end            
end