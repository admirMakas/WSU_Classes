K=[1 -1 0;
    -1 1 -1;
    0 -1 1];

x0 = [1;0;0];

for i = 1:10
   
    x_i = K*x0/norm(K*x0);
    
    lam = norm(K*x_i)
    
    x0 = x_i
    
end