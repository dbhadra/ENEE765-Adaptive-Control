clc;clear;
rt = linspace(0,1,100);
n = 100;

% r = sin(2*rt*2*pi);
% r = sin(rt*2*pi)+sin((2*rt+3)*2*pi)+sin((5*rt+4)*2*pi);
r = sin(2*rt*2*pi)+cos(rt*2*pi);

Rv1 = 0.6; 
Rv1_dB = 10*log10(Rv1); 
v = wgn(1,length(rt),Rv1_dB);

Rw1 = 0.3/n; 
Rw1_dB = 10*log10(Rw1); 
w = wgn(1,length(rt),Rw1_dB);


tspan = linspace(0,1,100);
mu0 = 0;
% opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t,mu] = ode45(@(t,mu) myode(t,mu,rt,r,v), tspan, mu0); 
plot(t,r,'k')
hold on;
plot(t,r+v,'b')
hold on;
plot(t,mu,'r')
xlabel('time (sec)','fontweight','bold','fontsize',16);
ylabel('test v/s measured v/s fused signal','fontweight','bold','fontsize',16);
legend('Test signal','Measured signal','Fused signal')

