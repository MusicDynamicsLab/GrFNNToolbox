function dB = Pa2dB(p)
    
    p_ref = 20/1000000;
    dB = 20*log10(p/p_ref);
