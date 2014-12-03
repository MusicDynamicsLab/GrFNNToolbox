function rt = relTimeHk(a, b, A, k, epsilon)
    
    rs = rstarHk(a, b, A, k, epsilon);
    a  = (a + 3*b*rs.^2);
    rt = -1./a;