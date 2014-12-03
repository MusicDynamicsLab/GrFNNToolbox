function ticks = tickMake(vector,numTicks)
len = length(vector);
things = linspace(eps('single'),1,numTicks);
ticks = vector(ceil(len*things));
if ticks(1) < 1
    coef = 10^(ceil(log10(ticks(1))) + 2);
    ticks = round(coef*ticks)/coef;
else
    ticks = round(ticks);
end
