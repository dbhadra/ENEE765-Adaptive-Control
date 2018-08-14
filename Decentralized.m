% Decentralized Adaptive Control of Interconnected Systems

clc
clear
k01 = 1.5; k02 = 1;
models = cell(1,2);

% Usual adaptive controller
models{1} = @(t,x) [(-2-x(5))*x(1)-(x(5)+1)*x(3)-0.5*(x(2)+x(4));... % Error e1
                     (-1.5-x(6))*x(2)-(x(6)+0.5)*x(4);... % Error e2
                     -x(3)+1;... % Local reference model xm1
                     -x(4)+4;... % Local reference model xm2
                     x(1)*(x(1)+x(3));... % Controller gain K1
                     x(2)*(x(2)+x(4))]; % Controller gain K2

% Modified controller
models{2} = @(t,x) [(-2-x(5))*x(1)-(x(5)+1)*x(3)-0.5*(x(2)+x(4));...
                       (-1.5-x(6))*x(2)-(x(6)+0.5)*x(4);...
                       -x(3)+1;...
                       -x(4)+4;...
                       x(1)*(x(1)+x(3))-0.05*x(5)*(abs(x(5)>k01));...
                       x(2)*(x(2)+x(4))-0.01*x(6)*(abs(x(6)>k02))];

for i=1:2
    
[t,xa] = ode45(models{i},[0 40],[0 0 0 0 0 0]);

figure(1)
subplot(2,1,2)
plot(t,xa(:,2))
xlabel('time','fontweight','bold','fontsize',16);
ylabel('error 2','fontweight','bold','fontsize',16);

subplot(2,1,1)
plot(t,xa(:,1))
xlabel('time','fontweight','bold','fontsize',16);
ylabel('error 1','fontweight','bold','fontsize',16);

figure(2)
subplot(2,1,1)
plot(t,xa(:,5))
xlabel('time','fontweight','bold','fontsize',16);
ylabel('K1(t)','fontweight','bold','fontsize',16);

subplot(2,1,2)
plot(t,xa(:,6))
xlabel('time','fontweight','bold','fontsize',16);
ylabel('K2(t)','fontweight','bold','fontsize',16);

figure(3)
subplot(2,1,1)
plot(t,xa(:,3))
title('xm1(t)')
xlabel('t'), ylabel('xm1')
subplot(2,1,2)
plot(t,xa(:,4))
title('xm2(t)')
xlabel('t'), ylabel('xm2')

pause(15)
end

% EOF