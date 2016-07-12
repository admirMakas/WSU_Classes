l=38
%define radii
r1 = 74/2
r2 = 61.21/2
r3 = 48.43/2
r4 = 35.64/2
r5 = 22.86/2
r6 = 10/2
%
%define areas
A1 = pi*r1^2
A2 = pi*r2^2
A3 = pi*r3^2
A4 = pi*r4^2
A5 = pi*r5^2
A6 = pi*r6^2
%
Ixx1 = (pi/4)*r1^4
Ixx2 = (pi/4)*r2^4
Ixx3 = (pi/4)*r3^4
Ixx4 = (pi/4)*r4^4
Ixx5 = (pi/4)*r5^4
Ixx6 = (pi/4)*r6^4
%
Iyy1 = Ixx1
Iyy2 = Ixx2
Iyy3 = Ixx3
Iyy4 = Ixx4
Iyy5 = Ixx5
Iyy6 = Ixx6
%
J1 = 0.95*(Ixx1+Iyy1)
J2 = 0.95*(Ixx2+Iyy2)
J3 = 0.95*(Ixx3+Iyy3)
J4 = 0.95*(Ixx4+Iyy4)
J5 = 0.95*(Ixx5+Iyy5)
J6 = 0.95*(Ixx6+Iyy6)