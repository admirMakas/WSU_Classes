function [Em] = getE(E, v)

    A=E/((1+v)*(1-2*v));
    
    Em=A*[1-v v v 0 0 0;
        v 1-v v 0 0 0;
        v v 1-v 0 0 0;
        0 0 0 (1-2*v)/2 0 0;
        0 0 0 0 (1-2*v)/2 0;
        0 0 0 0 0 (1-2*v)/2];

end