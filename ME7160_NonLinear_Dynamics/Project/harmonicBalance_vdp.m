%% <Define number of sample points
N = 100;
t = linspace(0,2*pi,N)';
t(end) = [];

%% <Defining function>
%x"+x'x^2-x'+x = cos(2t)sin(5t)
global func freq %Defining global variables
func = cos(t)+2*cos(5*t);

%% <Obtaining Frequency Spectrum>
dt = 1/N;
df=1./(N*dt); % Interval in freq
fpos=[0:N/2-1]'*df; % Freqs of positive half of spectrum
fneg=[N/2-1:-df:1]'*df*-1; %Freqs of negative half of spectrum
freq =[fpos;fneg]; %Full frequency spectrum

%% <Guessing intial solution for x>
x0 = ones(N-1,1)*2; %x=1

%% <Using Optimization toolbox to minimize residual and obtain x solution>
%See associated function 'resFUN_Duff.m'
options = optimset('MaxFunEvals',500000);
[x_HarBal,RES] = fminunc(@resFUN_vdp,x0,options);

%% < Using ODE45 to numerically solve for the ode
[t_num,x_num]=ode45('derivFUN_vdp',[0 100],[0 0]);

%% <Plotting solution with numerical solution>
t_shift = t+(1.935*pi);
figure
plot(t_num,x_num(:,1),t_shift,x_HarBal,'o')
axis([1.5*pi 4.5*pi -5 5]) %Shifting to steady state of solution
xlabel('t')
title('Harmoninc Balance -- Van Der Pol model')



