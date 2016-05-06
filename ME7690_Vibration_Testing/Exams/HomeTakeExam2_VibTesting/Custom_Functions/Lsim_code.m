clear all
clc

t=[0:0.001:10];
u=5*sin(10*t);

% define pertinent matrices for problem 2
M=[2 0; 0 7];
K=[2 -1; -1 2];
Cd=[0.1 -0.1; -0.1 0.2];

% define A matrix
MK=-K/M;
MC=-Cd/M;

display('A matrix below')
A=[zeros(2) eye(2); MK' MC']

% Define B matrix
Btilde=[1 -1];

B2= inv(M)*Btilde';

display('B matrix below')
B=[0 0 B2']'

C = [-1 0.5 -0.05 0.05]

D = 0.5

% system defined in state space, then it will be converted into a
% transfer function
sys = ss(A, B, C, D);

% define transfer function based on the state space model
H=tf(sys)

% define intial condition (not required to run the model)
xO = [0, 0, 0, 0];

y = lsim(H, u, t, xO);

plot(t, y)
hold on
plot(t, u)

%define nuerator and denominator of transfer function
num = [0.5, 0.007143, 0.07143, 0, 0];
den = [1, 0.07857, 1.286, 0.02857, 0.2143];

w1 = 0:0.001:6;
P = freqs(num, den, w1);
mag = abs(P);
phase = unwrap(angle(P));

phasedeg = phase*(180/pi);

figure
plot(w1, 20*log10(mag));
figure
plot(w1, phasedeg);