clear all;
clf;
format long;

a = 0;
b = 1;

%%%%% N = 2 %%%%%

N2 = 2;
h2 = (b - a)/N2;
Im2 = 0;

for i = 1:1:N2;
    xi = a + (i - 1/2)*h2;
    Im2 = Im2 + h2*exp(-xi^2);
end;

Im2

%%%%% N = 4 %%%%%

N4 = 4;
h4 = (b - a)/N4;
Im4 = 0;

for i = 1:1:N4;
    xi = a + (i - 1/2)*h4;
    Im4 = Im4 + h4*exp(-xi^2);
end;

Im4

%%%%% N = 8 %%%%%

N8 = 8;
h8 = (b - a)/N8;
Im8 = 0;

for i = 1:1:N8;
    xi = a + (i - 1/2)*h8;
    Im8 = Im8 + h8*exp(-xi^2);
end;

Im8

%%%%% N = 16 %%%%%

N16 = 16;
h16 = (b - a)/N16;
Im16 = 0;

for i = 1:1:N16;
    xi = a + (i - 1/2)*h16;
    Im16 = Im16 + h16*exp(-xi^2);
end;

Im16

%%%%% N = 32 %%%%%

N32 = 32;
h32 = (b - a)/N32;
Im32 = 0;

for i = 1:1:N32;
    xi = a + (i - 1/2)*h32;
    Im32 = Im32 + h32*exp(-xi^2);
end;

Im32

%%%%% N = 64 %%%%%

N64 = 64;
h64 = (b - a)/N64;
Im64 = 0;

for i = 1:1:N64;
    xi = a + (i - 1/2)*h64;
    Im64 = Im64 + h64*exp(-xi^2);
end;

Im64

%%%%% Richardson's Extrapolation %%%%%

IRE1 = 4/3*Im4 - 1/3*Im2
IRE2 = 4/3*Im8 - 1/3*Im4
IRE3 = 4/3*Im16 - 1/3*Im8
IRE4 = 4/3*Im32 - 1/3*Im16
IRE5 = 4/3*Im64 - 1/3*Im32

%%%%% Plotting %%%%%

Analytical = 1/2*sqrt(pi)*erf(b)
MidpointSubIntervalString = [N2,N4,N8,N16,N32,N64];
RESubIntervalString = [N4,N8,N16,N32,N64];
StepSizeString = [h2,h4,h8,h16,h32,h64];
MidpointErrorString = [Im2 - Analytical,Im4 - Analytical,Im8 - Analytical,Im16 - Analytical,Im32 - Analytical,Im64 - Analytical];
REErrorString = [IRE1 - Analytical,IRE2 - Analytical,IRE3 - Analytical,IRE4 - Analytical,IRE5 - Analytical];

loglog(MidpointSubIntervalString,MidpointErrorString,'-r');
hold on;
loglog(RESubIntervalString,-1*REErrorString,'-b');
loglog(MidpointSubIntervalString,0.1*StepSizeString.^(2),'-m');
loglog(MidpointSubIntervalString,0.025*(StepSizeString).^(4),'-g');
legend('Midpoint Error','Richardsons Extrapolation Error','h^2','h^4');
title('Error of Integration Approximations');
xlabel('Subintervals, N');
ylabel('Error');
hold off;