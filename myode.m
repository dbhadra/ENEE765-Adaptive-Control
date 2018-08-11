% Interpolation of a data for the sensor fusion consensus problem 

function dmudt = myode(t,mu,rt,r,v)
r = interp1(rt,r,t,'spline'); % Interpolate the data set (ft,f) at time t
v = interp1(rt,v,t,'spline'); % Interpolate the data set (gt,g) at time t
dmudt = 7*(r+v-mu);
end