t=0.01

A=[0 1; -10 0]

[V, D]=eig(A)

%next calculate matrix exponential of D
D=D*t
Dexp = expm(D)

Aexp=V*Dexp*inv(V)