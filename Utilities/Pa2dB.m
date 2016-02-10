%% Pa2dB
%  dB = Pa2dB(p)
%
%  Converts input in units of pascals p to output dB in dB SPL using the 
%  reference of 20 micropascals.

function dB = Pa2dB(p)
    
    p_ref = 20/1000000;
    dB = 20*log10(p/p_ref);
