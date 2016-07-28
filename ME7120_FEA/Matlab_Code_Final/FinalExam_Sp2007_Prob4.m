clear all
clc

%Generate reduced governing equations using Lagrange Mult. with 
%following constraint conditions.
%v1 = -v2 = (theta1*L)/(2^2)
%theta2 = theta1

%From the above conditions we have following C matrix that will be 
%appended to K and M in order to find reduced model.

syms theta E I m real
L=2;

C=[0.25 0 -0.25 theta];

K = ((E*I)/(L^3))*[12 6*L -12 6*L; 6*L 4*L^2 -6*L 2*L^2; -12 -6*L ...
    12 -6*L; 6*L 2*L^2 -6*L 4*L^2];

M = (m/420)*[156 22*L 54 -13*L; 22*L 4*L^2 13*L -3*L^2; 54 13*L 156 ...
    -22*L; -13*L -3*L^2 -22*L 4*L^2];

Kbar = [K C'];
Kbar = [Kbar; C, 0]

Mbar = [M C'];
Mbar = [Mbar; C, 0]

K11 = Kbar([2:2])
M11 = Mbar([2:2])

K12 = Kbar([1:1,3:5], [2:2])
M12 = Mbar([1:1,3:5], [2:2])

K21 = K12';
M21 = M12';

K22 = Kbar([1:1,3:5], [1:1,3:5])
M22 = Mbar([1:1,3:5], [1:1,3:5])

Q = [eye(1); -K22\K21]

