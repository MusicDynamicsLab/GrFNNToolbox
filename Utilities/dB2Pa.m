%% dB2Pa
%  p = dB2Pa(dB)
%
%  Converts input dB in dB SPL to units of pascals.

function p = dB2Pa(dB)
    
    p_ref = 20/1000000;
    p = p_ref*10.^(dB/20);
