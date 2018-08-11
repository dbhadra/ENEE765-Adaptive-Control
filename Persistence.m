% This contains different kinds of topologies with single/double integrator
% models with linear and nonlinear protocols

clc;clear;
k = 1;
D = [1 0 0 0 1;-1 1 0 -1 0;0 -1 1 0 -1;0 0 -1 1 0];
x0 = [0.1 0.1 0.77 0.8]';
y0 = [0.8 0.1 0.78 0]';
t = 0:0.1:5;
% g1 = (t.^2+1).*(t.^2+1).*(sin(t).^2);
L = zeros(4,4,length(t));
LapPotx = zeros(length(t),1);
LapPoty = zeros(length(t),1);


% for i=1:5
% g_sq(i,:) = (t.^2+1).*(t.^2+1).*(sin(i*t).^2);
% end

[t,x] = ode45(@(t,x) -k*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t.^2+1)^2*(sin(5*t)^2)])*D'*x,t, x0); 
[t,y] = ode45(@(t,y) -k*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2)])*D'*y,t, y0); 

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



% % Spanning tree edges
% 
xe = D'*x';
xe = xe';
xT = xe(:,1:3);
norm_xT = zeros(1,length(t));
for i =1:length(t)
    norm_xT(i) = norm(xe(i,:));
end
 plot(t,norm_xT') % Possible error in plot

% For n = 3 and triangle topology 

k = 1;
D = [-1 0 1;1 -1 0;0 1 -1];
x0 = [0.1 0.5 1]';
y0 = [1.5 0.9 2]';
t = 0:0.1:5;
LapPotx = zeros(length(t),1);
LapPoty = zeros(length(t),1);

% for i=1:5
% g_sq(i,:) = (t.^2+1).*(t.^2+1).*(sin(i*t).^2);
% end

[t,x] = ode45(@(t,x) -k*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D'*x,t, x0); 
[t,y] = ode45(@(t,y) -k*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D'*y,t, y0); 

T=0:0.1:5;
    for i =1:length(T),
      t=T(i);
L(:,:,i) = D*diag([(t.^2+1).*(t.^2+1).*(sin(t).^2),...
                   (t.^2+1).*(t.^2+1).*(sin(2*t).^2),...
                   (t.^2+1).*(t.^2+1).*(sin(3*t).^2)
              ])*D';
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

% For n = 6 and ~hexagon topology


k = 1;
D = [-1 0 0 0 -1 0 0 0 1;1 -1 0 0 0 0 0 0 -1;0 0 0 1 1 -1 0 0 0;0 1 -1 0 0 0 0 1 0;0 0 1 -1 0 1 -1 0 0;0 0 0 0 0 0 1 -1 0];
x0 = [0.1 0.1 1 2.5 0 2]';
y0 = [1.5 0.9 2 0 3 1]';
t = 0:0.1:5; % Start diverging with higher t? 

LapPotx = zeros(length(t),1);
LapPoty = zeros(length(t),1);
% for i=1:5
% g_sq(i,:) = (t.^2+1).*(t.^2+1).*(sin(i*t).^2);
% end

[t,x] = ode45(@(t,x) -D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D'*x,t, x0); 
[t,y] = ode45(@(t,y) -D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D'*y,t, y0); 

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

% For n = 3, triangle topology and double integrator dynamics 

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

% for i=1:5
% g_sq(i,:) = (t.^2+1).*(t.^2+1).*(sin(i*t).^2);
% end

[t,z] = ode45(@(t,z) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D' -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D'], 1)*z,t,z0); 
[t,w] = ode45(@(t,w) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D' -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2)])*D'], 1)*w,t,w0); 

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

% For n = 4, square topology and double integrator dynamics 
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
[t,z] = ode45(@(t,z) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2)])*D' -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2)])*D'], 1)*z,t,z0); 
[t,w] = ode45(@(t,w) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2)])*D' -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2)])*D'], 1)*w,t,w0); 
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

% Simple communication delay case
s = tf('s')
sys = exp(-0.1*s);    
sysx = pade(sys,3);
x0 = randn(10,1);
D = diag([2 3 4 4 4 4 4 4 3 2]);
A = [0 1 1 0 0 0 0 0 0 0;1 0 1 1 0 0 0 0 0 0;1 1 0 1 1 0 0 0 0 0;0 1 1 0 1 1 0 0 0 0;0 0 1 1 0 1 1 0 0 0;0 0 0 1 1 0 1 1 0 0;0 0 0 0 1 1 0 1 1 0;0 0 0 0 0 1 1 0 1 1;0 0 0 0 0 0 1 1 0 1;0 0 0 0 0 0 0 1 1 0];
L = D-A;
sol = dde23(@ddefun,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1],x0,[0 5]);
H = 1/(s+L*sysx);
K = H*x0;
[y,t] = impulse(K)


sol = dde23('ddefun',[1 1 1 1 1 1 1 1 1 1],[0 0 0 0 0 0 0 0 0 0],[0, 100]);

% Non-linear protocol case

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

% For n = 6 and ~hexagon topology, Double Integrator case

k = 1;
n = 6;
D = [-1 0 0 0 -1 0 0 0 1;1 -1 0 0 0 0 0 0 -1;0 0 0 1 1 -1 0 0 0;0 1 -1 0 0 0 0 1 0;0 0 1 -1 0 1 -1 0 0;0 0 0 0 0 0 1 -1 0];
x0 = [0.1 0.1 1 2.5 0 2]';
x0_dot = [0 0 0 0 0 0]';
y0 = [1.5 0.9 2 0 3 1]';
y0_dot = [0 0 0 0 0 0]';
z0 = [x0;x0_dot];
w0 = [y0;y0_dot];
t = 0:0.1:20; % Start diverging with higher t? 

g = 1;
[t,z] = ode45(@(t,z) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D' -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D'], 1)*z,t,z0); 
[t,w] = ode45(@(t,w) kron([zeros(n,n) eye(n);-D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D' -g*D*diag([(t^2+1)^2*(sin(t)^2) (t^2+1)^2*(sin(2*t)^2) (t^2+1)^2*(sin(3*t)^2) (t^2+1)^2*(sin(4*t)^2) (t^2+1)^2*(sin(5*t)^2) (t^2+1)^2*(sin(6*t)^2) (t^2+1)^2*(sin(7*t)^2) (t^2+1)^2*(sin(8*t)^2) (t^2+1)^2*(sin(9*t)^2)])*D'], 1)*w,t,w0);
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

% For n = 3 and triangle topology, constant weights

k = 1;
D = [-1 0 1;1 -1 0;0 1 -1];
x0 = [0.1 0.1 1]';
y0 = [1.5 0.9 2]';
t = 0:0.1:20;
W = eye(3);

% for i=1:5
% g_sq(i,:) = (t.^2+1).*(t.^2+1).*(sin(i*t).^2);
% end

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

% For n = 4 and square topology, constant weights

clc;clear;
k = 1;
D = [1 0 0 0 1;-1 1 0 -1 0;0 -1 1 0 -1;0 0 -1 1 0];
x0 = [0.1 0.1 0.77 0.78]';
y0 = [1.5 0.9 1.78 0]';
t = 0:0.1:20;
W = eye(5);

% for i=1:5
% g_sq(i,:) = (t.^2+1).*(t.^2+1).*(sin(i*t).^2);
% end

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

% Lyapunov analysis of spanning tree
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












