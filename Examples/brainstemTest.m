%% Basic network parameters
pmin = dB2Pa(  0); pmax = dB2Pa(120); 
e = 1/(pmax)^2;

sigFs   = 100000;

s = stimulusMake('fcn', [0 1], sigFs, {'exp'}, [300], dB2Pa(80), 0, ...
                 'ramp', 0.01, 1, 'display', 0);

brainstem;