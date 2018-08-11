
% Nonlinear protocol consensus problem

% function dxdt = vdp1(t,x)
% L = [2 -1 -1;0 1 -1;-1 0 1]; 
% dxdt = -L*[2*x(1)+sin(x(1)); 2*x(2)+sin(x(2)); 2*x(3)+sin(x(3))]; 
% end

function dxdt = vdp1(t,x)
L = [2 -1 -1;0 1 -1;-1 0 1]; 
dxdt = -L*[(x(1)>1)*x(1)^2+(x(1)>0 && x(1)<=1)*sqrt(x(1))+(x(1)>-1 && x(1)<=0)*sqrt(-x(1))+(x(1)<=-1)*-x(1)^2;...
    (x(2)>1)*x(2)^2+(x(2)>0 && x(2)<=1)*sqrt(x(2))+(x(2)>-1 && x(2)<=0)*sqrt(-x(2))+(x(2)<=-1)*-x(2)^2;...
        (x(3)>1)*x(3)^2+(x(3)>0 && x(3)<=1)*sqrt(x(3))+(x(3)>-1 && x(3)<=0)*sqrt(-x(3))+(x(3)<=-1)*-x(3)^2];
end