clc;clear;
D = [-1 0 1;0 1 -1;1 -1 0];
Le = D'*D;
u1 = 0:1:100;
k = 1;
Le_T = Le(1:2,1:2);
[V,D] = eig(Le_T);
u2 = 100;
eig_Le_T = eig(Le_T);
p = 2;
eta_sq = 2*k*min(eig_Le_T)*u1./((1+sqrt(p)*norm(D)*u2)^2);

a = ((p*u2)^1.5*max(eig_Le_T))./(1+2*u1*min(eig_Le_T));
b = 2*u1*min(eig_Le_T)./(max(eig_Le_T)*(1+2*u1*min(eig_Le_T)));
o = -a+sqrt(a.^2+b);
J = (1-o.^2/max(eig(inv(D))))./(1-eta_sq);

% For multiple values of u2
plot(u1,J)
axis([0 100 1 1.02])
hold on;

% Correct plot
xlabel('u1') % x-axis label
ylabel('J(p,u1,u2)') % y-axis label

