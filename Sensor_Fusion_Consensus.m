% Demonstrates how fusing data from multiple agents ultimately leads to correct consensus. Refer pages 6-7 in the PDF

clc;clear;
rt = linspace(0,1,100);
n = 100; % No. of measurements

r = sin(2*rt*2*pi)+cos(rt*2*pi); 
% Signal being measured by every sensor/node in a network/graph in a noisy environment.

Rv1 = 0.6; 
Rv1_dB = 10*log10(Rv1); 
v = wgn(1,length(rt),Rv1_dB); % Generates white Gaussian noise of Rv1_dB power

Rw1 = 0.3/n; 
Rw1_dB = 10*log10(Rw1); 
w = wgn(1,length(rt),Rw1_dB);

tspan = linspace(0,1,100);
mu0 = 0;
[t,mu] = ode45(@(t,mu) myode(t,mu,rt,r,v), tspan, mu0); % Refer myode.m
plot(t,r,'k')
hold on;
plot(t,r+v,'b')
plot(t,mu,'r')
xlabel('time (sec)','fontweight','bold','fontsize',16);
ylabel('test v/s measured v/s fused signal','fontweight','bold','fontsize',16);
legend('Test signal','Measured signal','Fused signal')
hold off;

% EOF
