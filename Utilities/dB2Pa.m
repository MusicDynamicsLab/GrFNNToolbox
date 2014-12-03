function p = dB2Pa(dB)
    
    p_ref = 20/1000000;
    p = p_ref*10.^(dB/20);
