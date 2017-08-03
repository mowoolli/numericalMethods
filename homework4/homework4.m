clear all;
format long;

h1 = sqrt(1/3480);
h2 = sqrt(1/880);
r = h2/h1;

load Gauss20x20;
load Gauss40x40;
load leastSquares20x20;
load leastSquares40x40;
load Gauss20x20uniform;
load Gauss40x40uniform;
load leastSquares20x20uniform;
load leastSquares40x40uniform;

y = 0:0.0000001:1;
Analytical = sin(pi()/2*y);
analyticalGradient = pi()/2*cos(pi()/2*y);

%%%%% PLOTTING %%%%%

% plot(Gauss20x20(:,1),Gauss20x20(:,2),'*r',y,Analytical,'-g');
% set(gca,'FontSize',16);
% xlabel('y','FontSize',16);
% ylabel('T(y)','FontSize',16);
% title('Solution on 20x20','FontSize',16);
% legend('OpenFOAM','Analytical');
% pause;
% 
% plot(Gauss20x20(:,1),Gauss20x20(:,3),'*r',leastSquares20x20(:,1),leastSquares20x20(:,3),'*b',y,analyticalGradient,'-g');
% set(gca,'FontSize',16);
% xlabel('y','FontSize',16);
% ylabel('Grad(T)','FontSize',16);
% title('Solution of Gradient on 20x20','FontSize',16);
% legend('Gauss','Least Squares','Analytical');
% pause;
% 
% plot(Gauss40x40(:,1),Gauss40x40(:,2),'*r',y,Analytical,'-g');
% set(gca,'FontSize',16);
% xlabel('y','FontSize',16);
% ylabel('T(y)','FontSize',16);
% title('Solution on 40x40','FontSize',16);
% legend('OpenFOAM','Analytical');
% pause;
% 
% plot(Gauss40x40(:,1),Gauss40x40(:,3),'*r',leastSquares40x40(:,1),leastSquares40x40(:,3),'*b',y,analyticalGradient,'-g');
% set(gca,'FontSize',16);
% xlabel('y','FontSize',16);
% ylabel('Grad(T)','FontSize',16);
% title('Solution of Gradient on 40x40','FontSize',16);
% legend('Gauss','Least Squares','Analytical');
% pause;
% 
% plot(Gauss20x20uniform(:,1),Gauss20x20uniform(:,2),'*r',y,Analytical,'-g');
% set(gca,'FontSize',16);
% xlabel('y','FontSize',16);
% ylabel('T(y)','FontSize',16);
% title('Solution on 20x20','FontSize',16);
% legend('OpenFOAM','Analytical');
% pause;
% 
% plot(Gauss20x20uniform(:,1),Gauss20x20uniform(:,3),'*r',leastSquares20x20uniform(:,1),leastSquares20x20uniform(:,3),'*b',y,analyticalGradient,'-g');
% set(gca,'FontSize',16);
% xlabel('y','FontSize',16);
% ylabel('Grad(T)','FontSize',16);
% title('Solution of Gradient on 20x20','FontSize',16);
% legend('Gauss','Least Squares','Analytical');
% pause;
% 
% plot(Gauss40x40uniform(:,1),Gauss40x40uniform(:,2),'*r',y,Analytical,'-g');
% set(gca,'FontSize',16);
% xlabel('y','FontSize',16);
% ylabel('T(y)','FontSize',16);
% title('Solution on 40x40','FontSize',16);
% legend('OpenFOAM','Analytical');
% pause;
% 
% plot(Gauss40x40uniform(:,1),Gauss40x40uniform(:,3),'*r',leastSquares40x40uniform(:,1),leastSquares40x40uniform(:,3),'*b',y,analyticalGradient,'-g');
% set(gca,'FontSize',16);
% xlabel('y','FontSize',16);
% ylabel('Grad(T)','FontSize',16);
% title('Solution of Gradient on 40x40','FontSize',16);
% legend('Gauss','Least Squares','Analytical');
% pause;

%%%%% OBSERVED ORDER OF ACCURACY %%%%%

errorGauss20x20 = zeros(length(Gauss20x20uniform(:,1)),1);
errorLeastSquares20x20 = zeros(length(Gauss20x20uniform(:,1)),1);
errorGauss40x40 = zeros(length(Gauss40x40uniform(:,1)),1);
errorLeastSquares40x40 = zeros(length(Gauss40x40uniform(:,1)),1);

for i = 1:length(Gauss20x20uniform(:,1));
    errorGauss20x20(i) = Gauss20x20uniform(i,3) - pi()/2*cos(pi()/2*Gauss20x20uniform(i,1));
    errorLeastSquares20x20(i) = leastSquares20x20uniform(i,3) - pi()/2*cos(pi()/2*leastSquares20x20uniform(i,1));
    errorGauss40x40(i) = Gauss40x40uniform(i,3) - pi()/2*cos(pi()/2*Gauss40x40uniform(i,1));
    errorLeastSquares40x40(i) = leastSquares40x40uniform(i,3) - pi()/2*cos(pi()/2*leastSquares40x40uniform(i,1));
end

observedOrderGauss = mean(abs(log(abs(errorGauss20x20./errorGauss40x40))./log(r)));
observedOrderLeastSquares = mean(abs(log(abs(errorLeastSquares20x20./errorLeastSquares40x40))./log(r)));



extrapolatedGauss = (r^2/(r^2-1))*Gauss20x20uniform - Gauss40x40uniform/(r^2-1);
extrapolatedLeastSquares = (r^2/(r^2-1))*leastSquares20x20uniform - leastSquares40x40uniform/(r^2-1);

approximateRelativeErrorGauss = abs((Gauss20x20uniform(:,3) - Gauss40x40uniform(:,3))./Gauss20x20uniform(:,3));
approximateRelativeErrorLeastSquares = abs((leastSquares20x20uniform(:,3) - leastSquares40x40uniform(:,3))./leastSquares20x20uniform(:,3));

GCIGauss = 1.25.*approximateRelativeErrorGauss/(r^(observedOrderGauss) - 1);
GCILeastSquares = 1.25.*approximateRelativeErrorLeastSquares/(r^(observedOrderLeastSquares) - 1);

errorbar(Gauss40x40uniform(:,1),Gauss40x40uniform(:,3),GCIGauss,'*r');
hold on;
errorbar(leastSquares40x40uniform(:,1),leastSquares40x40uniform(:,3),GCILeastSquares,'*b');
plot(y,analyticalGradient,'-g');
set(gca,'FontSize',16);
xlabel('y','FontSize',16);
ylabel('Grad(T)','FontSize',16);
title('Solution of Gradient on 40x40','FontSize',16);
legend('Gauss','Least Squares','Analytical');
axis([0 1 0 1.6]);
hold off;
