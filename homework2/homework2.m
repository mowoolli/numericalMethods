clear all;
format long;
clf;

a = 0;
b = 1;

%%%%% N = 10 %%%%%

N = 10;
h = (b - a)/N;
u10Euler(1) = 1;
u10Heuns(1) = 1;
u10RK4(1) = 1;
u10AB2(1) = 1;

for i = 1:N;
        u10Euler(i+1) = u10Euler(i) + h*(-u10Euler(i)*tan(h*(i-1)) + sin(2*h*(i-1)));
        u10Heuns(i+1) = u10Heuns(i) + (h/2)*((-u10Heuns(i)*tan(h*(i-1)) + sin(2*h*(i-1))) + (-(u10Heuns(i) + h*(-u10Heuns(i)*tan(h*(i-1)) + sin(2*h*(i-1))))*tan(h*(i)) + sin(2*h*(i))));
        k1(i) = h*(-u10RK4(i)*tan(h*(i-1)) + sin(2*h*(i-1)));
        k2(i) = h*(-(u10RK4(i) + 0.5*k1(i))*tan(h*(i-1)+h/2) + sin(2*(h*(i-1)+h/2)));
        k3(i) = h*(-(u10RK4(i) + 0.5*k2(i))*tan(h*(i-1)+h/2) + sin(2*(h*(i-1) + h/2)));
        k4(i) = h*(-(u10RK4(i) + k3(i))*tan(h*(i-1)+h) + sin(2*(h*(i-1) + h)));
        u10RK4(i+1) = u10RK4(i) + 1/6*(k1(i) + 2*k2(i) + 2*k3(i) + k4(i));
        u10AB2(2) = u10RK4(2);
        u10AB2(i+2) = u10AB2(i+1) + (h/2)*(3*(-u10AB2(i+1)*tan(h*(i)) + sin(2*h*(i)))-(-u10AB2(i)*tan(h*(i-1)) + sin(2*h*(i-1))));
end

%%%%% N = 20 %%%%%

N = 20;
h = (b - a)/N;
u20Euler(1) = 1;
u20Heuns(1) = 1;
u20RK4(1) = 1;
u20AB2(1) = 1;

for i = 1:N;
        u20Euler(i+1) = u20Euler(i) + h*(-u20Euler(i)*tan(h*(i-1)) + sin(2*h*(i-1)));
        u20Heuns(i+1) = u20Heuns(i) + (h/2)*((-u20Heuns(i)*tan(h*(i-1)) + sin(2*h*(i-1))) + (-(u20Heuns(i) + h*(-u20Heuns(i)*tan(h*(i-1)) + sin(2*h*(i-1))))*tan(h*(i)) + sin(2*h*(i)))); 
        k1(i) = h*(-u20RK4(i)*tan(h*(i-1)) + sin(2*h*(i-1)));
        k2(i) = h*(-(u20RK4(i) + 0.5*k1(i))*tan(h*(i-1)+h/2) + sin(2*(h*(i-1)+h/2)));
        k3(i) = h*(-(u20RK4(i) + 0.5*k2(i))*tan(h*(i-1)+h/2) + sin(2*(h*(i-1) + h/2)));
        k4(i) = h*(-(u20RK4(i) + k3(i))*tan(h*(i-1)+h) + sin(2*(h*(i-1) + h)));
        u20RK4(i+1) = u20RK4(i) + 1/6*(k1(i) + 2*k2(i) + 2*k3(i) + k4(i));
        u20AB2(2) = u20RK4(2);
        u20AB2(i+2) = u20AB2(i+1) + (h/2)*(3*(-u20AB2(i+1)*tan(h*(i)) + sin(2*h*(i)))-(-u20AB2(i)*tan(h*(i-1)) + sin(2*h*(i-1))));
end

%%%%% N = 40 %%%%%

N = 40;
h = (b - a)/N;
u40Euler(1) = 1;
u40Heuns(1) = 1;
u40RK4(1) = 1;
u40RK4(1) = 1;
u40AB2(1) = 1;

for i = 1:N;
        u40Euler(i+1) = u40Euler(i) + h*(-u40Euler(i)*tan(h*(i-1)) + sin(2*h*(i-1)));
        u40Heuns(i+1) = u40Heuns(i) + (h/2)*((-u40Heuns(i)*tan(h*(i-1)) + sin(2*h*(i-1))) + (-(u40Heuns(i) + h*(-u40Heuns(i)*tan(h*(i-1)) + sin(2*h*(i-1))))*tan(h*(i)) + sin(2*h*(i))));  
        k1(i) = h*(-u40RK4(i)*tan(h*(i-1)) + sin(2*h*(i-1)));
        k2(i) = h*(-(u40RK4(i) + 0.5*k1(i))*tan(h*(i-1)+h/2) + sin(2*(h*(i-1)+h/2)));
        k3(i) = h*(-(u40RK4(i) + 0.5*k2(i))*tan(h*(i-1)+h/2) + sin(2*(h*(i-1) + h/2)));
        k4(i) = h*(-(u40RK4(i) + k3(i))*tan(h*(i-1)+h) + sin(2*(h*(i-1) + h)));
        u40RK4(i+1) = u40RK4(i) + 1/6*(k1(i) + 2*k2(i) + 2*k3(i) + k4(i));
        u40AB2(2) = u40RK4(2);
        u40AB2(i+2) = u40AB2(i+1) + (h/2)*(3*(-u40AB2(i+1)*tan(h*(i)) + sin(2*h*(i)))-(-u40AB2(i)*tan(h*(i-1)) + sin(2*h*(i-1))));
end

%%%%% N = 80 %%%%%

N = 80;
h = (b - a)/N;
u80Euler(1) = 1;
u80Heuns(1) = 1;
u80RK4(1) = 1;
u80RK4(1) = 1;
u80AB2(1) = 1;

for i = 1:N;
        u80Euler(i+1) = u80Euler(i) + h*(-u80Euler(i)*tan(h*(i-1)) + sin(2*h*(i-1)));
        u80Heuns(i+1) = u80Heuns(i) + (h/2)*((-u80Heuns(i)*tan(h*(i-1)) + sin(2*h*(i-1))) + (-(u80Heuns(i) + h*(-u80Heuns(i)*tan(h*(i-1)) + sin(2*h*(i-1))))*tan(h*(i)) + sin(2*h*(i))));    
        k1(i) = h*(-u80RK4(i)*tan(h*(i-1)) + sin(2*h*(i-1)));
        k2(i) = h*(-(u80RK4(i) + 0.5*k1(i))*tan(h*(i-1)+h/2) + sin(2*(h*(i-1)+h/2)));
        k3(i) = h*(-(u80RK4(i) + 0.5*k2(i))*tan(h*(i-1)+h/2) + sin(2*(h*(i-1) + h/2)));
        k4(i) = h*(-(u80RK4(i) + k3(i))*tan(h*(i-1)+h) + sin(2*(h*(i-1) + h)));
        u80RK4(i+1) = u80RK4(i) + 1/6*(k1(i) + 2*k2(i) + 2*k3(i) + k4(i));
        u80AB2(2) = u80RK4(2);
        u80AB2(i+2) = u80AB2(i+1) + (h/2)*(3*(-u80AB2(i+1)*tan(h*(i)) + sin(2*h*(i)))-(-u80AB2(i)*tan(h*(i-1)) + sin(2*h*(i-1))));
end

%%%%% N = 160 %%%%%

N = 160;
h = (b - a)/N;
u160Euler(1) = 1;
u160Heuns(1) = 1;
u160RK4(1) = 1;
u160RK4(1) = 1;
u160AB2(1) = 1;

for i = 1:N;
        u160Euler(i+1) = u160Euler(i) + h*(-u160Euler(i)*tan(h*(i-1)) + sin(2*h*(i-1)));
        u160Heuns(i+1) = u160Heuns(i) + (h/2)*((-u160Heuns(i)*tan(h*(i-1)) + sin(2*h*(i-1))) + (-(u160Heuns(i) + h*(-u160Heuns(i)*tan(h*(i-1)) + sin(2*h*(i-1))))*tan(h*(i)) + sin(2*h*(i))));       
        k1(i) = h*(-u160RK4(i)*tan(h*(i-1)) + sin(2*h*(i-1)));
        k2(i) = h*(-(u160RK4(i) + 0.5*k1(i))*tan(h*(i-1)+h/2) + sin(2*(h*(i-1)+h/2)));
        k3(i) = h*(-(u160RK4(i) + 0.5*k2(i))*tan(h*(i-1)+h/2) + sin(2*(h*(i-1) + h/2)));
        k4(i) = h*(-(u160RK4(i) + k3(i))*tan(h*(i-1)+h) + sin(2*(h*(i-1) + h)));
        u160RK4(i+1) = u160RK4(i) + 1/6*(k1(i) + 2*k2(i) + 2*k3(i) + k4(i));
        u160AB2(2) = u160RK4(2);
        u160AB2(i+2) = u160AB2(i+1) + (h/2)*(3*(-u160AB2(i+1)*tan(h*(i)) + sin(2*h*(i)))-(-u160AB2(i)*tan(h*(i-1)) + sin(2*h*(i-1))));
end

%%%%% PLOTTING %%%%%

u10AB2 = u10AB2(1:11);
u20AB2 = u20AB2(1:21);
u40AB2 = u40AB2(1:41);
u80AB2 = u80AB2(1:81);
u160AB2 = u160AB2(1:161);

t = 0:0.00001:1;
Analytical = 3*cos(t) - 2*cos(t).^2;

EulerErrorString = [u10Euler(end)-Analytical(end) u20Euler(end)-Analytical(end) u40Euler(end)-Analytical(end) u80Euler(end)-Analytical(end) u160Euler(end)-Analytical(end)];
HeunsErrorString = [u10Heuns(end)-Analytical(end) u20Heuns(end)-Analytical(end) u40Heuns(end)-Analytical(end) u80Heuns(end)-Analytical(end) u160Heuns(end)-Analytical(end)];
RK4ErrorString = [u10RK4(end)-Analytical(end) u20RK4(end)-Analytical(end) u40RK4(end)-Analytical(end) u80RK4(end)-Analytical(end) u160RK4(end)-Analytical(end)];
AB2ErrorString = [u10AB2(end)-Analytical(end) u20AB2(end)-Analytical(end) u40AB2(end)-Analytical(end) u80AB2(end)-Analytical(end) u160AB2(end)-Analytical(end)];
SubintervalString = [10 20 40 80 160];
StepSizeString = [(b-a)/10 (b-a)/20 (b-a)/40 (b-a)/80 (b-a)/160];
time = 0:(b-a)/10:1;

plot(time,u10Euler,'*m');
hold on;
plot(time,u10Heuns,'*r');
plot(time,u10RK4,'*b');
plot(time,u10AB2,'*g');
plot(t,Analytical,'-k');
legend('Euler','Heuns','RK4','AB2','Analytical');
title('Comparison of Approximations and Analytical');
xlabel('Time, t');
ylabel('Velocity, u');
hold off;
pause

time = 0:(b-a)/40:1;

plot(time,u40Euler,time,u40Heuns,time,u40RK4,time,u40AB2,t,Analytical);
pause

time = 0:(b-a)/160:1;
plot(time,u160Euler,time,u160Heuns,time,u160RK4,time,u160AB2,t,Analytical);
pause
loglog(SubintervalString,abs(EulerErrorString),'*m');
hold on;
loglog(SubintervalString,abs(HeunsErrorString),'*r');
loglog(SubintervalString,abs(RK4ErrorString),'*b');
loglog(SubintervalString,abs(AB2ErrorString),'*g');
loglog(SubintervalString,StepSizeString.^1,'--k');
loglog(SubintervalString,0.1*StepSizeString.^2,'-.k');
loglog(SubintervalString,0.8*StepSizeString.^4,'-k');
legend('Euler','Heuns','RK4','AB2','h','h^2','h^4');
title('Error of Approximations');
xlabel('Subintervals, N');
ylabel('|Error|');
hold off;