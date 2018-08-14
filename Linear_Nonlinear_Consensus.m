clc;clear;

%% Consensus with Linear Protocols

k = 1; % control gain
D = [1 0 0 0 1;...
    -1 1 0 -1 0;...
    0 -1 1 0 -1;...
    0 0 -1 1 0]; % incidence matrix
x0 = [0.1 0.1 0.77 0.8]';
y0 = [0.8 0.1 0.78 0]';
t = 0:0.1:5;
L = zeros(4,4,length(t));
LapPotx = zeros(length(t),1);
LapPoty = zeros(length(t),1);
% Laplacian potential as a measure of total disagreement among the nodes

[t,x] = ode45(@(t,x) -k*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2)... 
                (t.^2+1)^2*(sin(5*t)^2)])*D'*x,t, x0); 
            
[t,y] = ode45(@(t,y) -k*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2)...
                (t^2+1)^2*(sin(5*t)^2)])*D'*y,t, y0); 
            
% (x,y) State trajectory of the 4 agents in the network. Refer (9) in PDF          

T=0:0.1:5;
for i =1:length(T),
t=T(i);
L(:,:,i) = D*diag([(t.^2+1).*(t.^2+1).*(sin(t).^2),...
                (t.^2+1).*(t.^2+1).*(sin(2*t).^2),...
                (t.^2+1).*(t.^2+1).*(sin(3*t).^2),...
                (t.^2+1).*(t.^2+1).*(sin(4*t).^2),...
                (t.^2+1).*(t.^2+1).*(sin(5*t).^2)])*D';
            
LapPotx(i,1) = x(i,:)*L(:,:,i)*x(i,:)';
LapPoty(i,1) = y(i,:)*L(:,:,i)*y(i,:)';

end
% L is the Laplacian Matrix
   
    
subplot(2,1,1)
plot(T,LapPotx)
xlabel('time','fontweight','bold','fontsize',16);
title('Laplacian Potential of x-coordinates','fontweight','bold','fontsize',16);
subplot(2,1,2)
plot(T,LapPoty)
xlabel('time','fontweight','bold','fontsize',16);
title('Laplacian Potential of y-coordinates','fontweight','bold','fontsize',16);

figure(1)
plot(x(:,1),y(:,1),'b')
hold on;
plot(x(:,2),y(:,2),'g')
hold on;
plot(x(:,3),y(:,3),'r')
hold on;
plot(x(:,4),y(:,4),'k')
xlabel('x-coordinates','fontweight','bold','fontsize',16);
ylabel('y-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3','Agent 4')

pause(15)

%% For n = 3 and triangle topology (n : no. of agents)
clc;clear;
k = 1;
D = [-1 0 1;1 -1 0;0 1 -1];
x0 = [0.1 0.5 1]';
y0 = [1.5 0.9 2]';
t = 0:0.1:5;
LapPotx = zeros(length(t),1);
LapPoty = zeros(length(t),1);

[t,x] = ode45(@(t,x) -k*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D'*x,t, x0); 
[t,y] = ode45(@(t,y) -k*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D'*y,t, y0); 

T=0:0.1:5;
for i =1:length(T),
t=T(i);
L(:,:,i) = D*diag([(t.^2+1).*(t.^2+1).*(sin(t).^2),(t.^2+1).*(t.^2+1).*(sin(2*t).^2),(t.^2+1).*(t.^2+1).*(sin(3*t).^2)])*D';
LapPotx(i,1) = x(i,:)*L(:,:,i)*x(i,:)';
LapPoty(i,1) = y(i,:)*L(:,:,i)*y(i,:)';
end
   
    
subplot(2,1,1)
plot(T,LapPotx)
xlabel('time','fontweight','bold','fontsize',16);
title('Laplacian Potential of x-coordinates','fontweight','bold','fontsize',16);
subplot(2,1,2)
plot(T,LapPoty)
xlabel('time','fontweight','bold','fontsize',16);
title('Laplacian Potential of y-coordinates','fontweight','bold','fontsize',16);

figure(1)
plot(x(:,1),y(:,1),'k')
hold on;
plot(x(:,2),y(:,2),'g')
hold on;
plot(x(:,3),y(:,3),'r')
xlabel('x-coordinates','fontweight','bold','fontsize',16);
ylabel('y-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3')

pause(15)

%% For n = 6 and hexagon topology
clc;clear;

k = 1;
D = [-1 0 0 0 -1 0 0 0 1;...
    1 -1 0 0 0 0 0 0 -1;...
    0 0 0 1 1 -1 0 0 0;...
    0 1 -1 0 0 0 0 1 0;...
    0 0 1 -1 0 1 -1 0 0;...
    0 0 0 0 0 0 1 -1 0];
x0 = [0.1 0.1 1 2.5 0 2]';
y0 = [1.5 0.9 2 0 3 1]';
t = 0:0.1:5; 

LapPotx = zeros(length(t),1);
LapPoty = zeros(length(t),1);

[t,x] = ode45(@(t,x) -D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2)...
                (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D'*x,t, x0); 
            
[t,y] = ode45(@(t,y) -D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2)...
                (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D'*y,t, y0); 

T=0:0.1:5;
for i =1:length(T),
t=T(i);
L(:,:,i) = D*diag([(t.^2+1).*(t.^2+1).*(sin(t).^2),...
                   (t.^2+1).*(t.^2+1).*(sin(2*t).^2),...
                   (t.^2+1).*(t.^2+1).*(sin(3*t).^2),...
                   (t.^2+1).*(t.^2+1).*(sin(4*t).^2),... 
                   (t.^2+1).*(t.^2+1).*(sin(5*t).^2),...
                   (t.^2+1).*(t.^2+1).*(sin(6*t).^2),...
                   (t.^2+1).*(t.^2+1).*(sin(7*t).^2),...
                   (t.^2+1).*(t.^2+1).*(sin(8*t).^2),...
                   (t.^2+1).*(t.^2+1).*(sin(9*t).^2)])*D';
LapPotx(i,1) = x(i,:)*L(:,:,i)*x(i,:)';
LapPoty(i,1) = y(i,:)*L(:,:,i)*y(i,:)';

end
   
    
subplot(2,1,1)
plot(T,LapPotx)
xlabel('time','fontweight','bold','fontsize',16);
title('Laplacian Potential of x-coordinates','fontweight','bold','fontsize',16);
subplot(2,1,2)
plot(T,LapPoty)
xlabel('time','fontweight','bold','fontsize',16);
title('Laplacian Potential of y-coordinates','fontweight','bold','fontsize',16);


figure(1)
plot(x(:,1),y(:,1),'c')
hold on;
plot(x(:,2),y(:,2),'g')
hold on;
plot(x(:,3),y(:,3),'r')
hold on;
plot(x(:,4),y(:,4) ,'b')
hold on;
plot(x(:,5),y(:,5),'m')
hold on;
plot(x(:,6),y(:,6),'k')
xlabel('x-coordinates','fontweight','bold','fontsize',16);
ylabel('y-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6')

pause(15)

%% For n = 3, triangle topology and double integrator dynamics 
clc;clear;

k = 1;
n = 3;
D = [-1 0 1;1 -1 0;0 1 -1];
x0 = [0.1 0.1 1]';
x0_dot = [0 0 0]';
y0 = [1.5 0.9 2]';
y0_dot = [0 0 0]';
z0 = [x0;x0_dot];
w0 = [y0;y0_dot];
t = 0:0.1:20;
g=1;

[t,z] = ode45(@(t,z) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D'...
                -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D'], 1)*z,t,z0); 
            
[t,w] = ode45(@(t,w) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D'...
                -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D'], 1)*w,t,w0); 

% kron : Kronecker delta function            
            
x = z(:,1:3);
y = w(:,1:3);
figure(1)
plot(x(:,1),y(:,1),'k')
hold on;
plot(x(:,2),y(:,2),'g')
hold on;
plot(x(:,3),y(:,3),'r')
xlabel('x-coordinates','fontweight','bold','fontsize',16);
ylabel('y-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3')

pause(15)

%% For n = 4, square topology and double integrator dynamics 
clc;clear;

k = 1;
n = 4;
D = [1 0 0 0 1;-1 1 0 -1 0;0 -1 1 0 -1;0 0 -1 1 0];
x0 = [0.1 0.1 0.77 0.78]';
x0_dot = [0 0 0 0]';
y0 = [1.5 0.9 1.78 0]';
y0_dot = [0 0 0 0]';
z0 = [x0;x0_dot];
w0 = [y0;y0_dot];
t = 0:0.1:20;
g = 1;

[t,z] = ode45(@(t,z) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2)])*D'...
                -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2)])*D'], 1)*z,t,z0); 
            
[t,w] = ode45(@(t,w) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2)])*D'...
                -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2)])*D'], 1)*w,t,w0);
            
x = z(:,1:4);
y = w(:,1:4);
figure(1)
plot(x(:,1),y(:,1),'b')
hold on;
plot(x(:,2),y(:,2),'g')
hold on;
plot(x(:,3),y(:,3),'r')
hold on;
plot(x(:,4),y(:,4),'k')
xlabel('x-coordinates','fontweight','bold','fontsize',16);
ylabel('y-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3','Agent 4')

pause(15)

%% Consensus with nonlinear protocol
clc;clear;

x0 = [0.1 0.5 1]';
y0 = [1.5 0.9 2]';
t = 0:0.01:5;
[t,x] = ode45(@vdp1,t,x0);
[t,y] = ode45(@vdp1,t,y0);

figure(1)
plot(x(:,1),y(:,1),'k')
hold on;
plot(x(:,2),y(:,2),'g')
hold on;
plot(x(:,3),y(:,3),'r')
xlabel('x-coordinates','fontweight','bold','fontsize',16);
ylabel('y-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3')
subplot(2,1,1)
plot(t,x)
xlabel('time','fontweight','bold','fontsize',16);
ylabel('x-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3')
subplot(2,1,2)
plot(t,y)
xlabel('time','fontweight','bold','fontsize',16);
ylabel('y-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3')

pause(15)

%% For n = 6 and ~hexagon topology, Double Integrator case
clc;clear;

k = 1;
n = 6;
D = [-1 0 0 0 -1 0 0 0 1;...
    1 -1 0 0 0 0 0 0 -1;...
    0 0 0 1 1 -1 0 0 0;...
    0 1 -1 0 0 0 0 1 0;...
    0 0 1 -1 0 1 -1 0 0;...
    0 0 0 0 0 0 1 -1 0];
x0 = [0.1 0.1 1 2.5 0 2]';
x0_dot = [0 0 0 0 0 0]';
y0 = [1.5 0.9 2 0 3 1]';
y0_dot = [0 0 0 0 0 0]';
z0 = [x0;x0_dot];
w0 = [y0;y0_dot];
t = 0:0.1:20; 

g = 1;
[t,z] = ode45(@(t,z) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)...
                (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D' -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D'], 1)*z,t,z0); 

[t,w] = ode45(@(t,w) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)...
                (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D' -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D'], 1)*w,t,w0);

figure(1)
plot(x(:,1),y(:,1),'c')
hold on;
plot(x(:,2),y(:,2),'g')
hold on;
plot(x(:,3),y(:,3),'r')
hold on;
plot(x(:,4),y(:,4) ,'b')
hold on;
plot(x(:,5),y(:,5),'m')
hold on;
plot(x(:,6),y(:,6),'k')
xlabel('x-coordinates','fontweight','bold','fontsize',16);
ylabel('y-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6')

pause(15)

%% For n = 3 and triangle topology, constant weights
clc;clear;

k = 1;
D = [-1 0 1;1 -1 0;0 1 -1];
x0 = [0.1 0.1 1]';
y0 = [1.5 0.9 2]';
t = 0:0.1:20;
W = eye(3);

[t,x] = ode45(@(t,x) -D*W*D'*x,t, x0); 
[t,y] = ode45(@(t,y) -D*W*D'*y,t, y0); 
figure(1)
plot(x(:,1),y(:,1),'k')
hold on;
plot(x(:,2),y(:,2),'g')
hold on;
plot(x(:,3),y(:,3),'r')
xlabel('x-coordinates','fontweight','bold','fontsize',16);
ylabel('y-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3')

pause(15)

%% For n = 4 and square topology, constant weights
clc;clear;

k = 1;
D = [1 0 0 0 1;-1 1 0 -1 0;0 -1 1 0 -1;0 0 -1 1 0];
x0 = [0.1 0.1 0.77 0.78]';
y0 = [1.5 0.9 1.78 0]';
t = 0:0.1:20;
W = eye(5);

[t,x] = ode45(@(t,x) -k*D*W*D'*x,t, x0); 
[t,y] = ode45(@(t,y) -k*D*W*D'*y,t, y0); 
figure(1)
plot(x(:,1),y(:,1),'b')
hold on;
plot(x(:,2),y(:,2),'g')
hold on;
plot(x(:,3),y(:,3),'r')
hold on;
plot(x(:,4),y(:,4),'k')
xlabel('x-coordinates','fontweight','bold','fontsize',16);
ylabel('y-coordinates','fontweight','bold','fontsize',16);
legend('Agent 1','Agent 2','Agent 3','Agent 4')

pause(15)

%% Lyapunov analysis of spanning tree
clc;clear;

k = 1;
DT = [1 0 0;1 1 0;0 1 1;0 0 1];
x0 = [0.1 0.5 0.9]';
y0 = [0.9 0.5 0.1]';
t = 0:0.1:20;
LT = DT'*DT;
eigLT = inv(diag([eig(LT)])); 
eigV = zeros(length(t),3);

for i =1:length(t)
   v = eigLT*LT*diag([(t.^2+1).*(t.^2+1).*(sin(t).^2) (t.^2+1).*(t.^2+1).*(sin(2*t).^2) (t.^2+1).*(t.^2+1).*(sin(3*t).^2)]);
   eigV(i,:) = eig(v+v'); 
end

% EOF










