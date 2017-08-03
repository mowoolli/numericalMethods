clear all;
format long;

load Gauss20x20neumann;
load LeastSquares20x20neumann;
load Gauss40x40neumann;
load LeastSquares40x40neumann;

a = 1;
b = 1;

Nx = 61;
Ny = 61;

deltax = a/(Nx-1);
deltay = b/(Ny-1);

K = 101;

uGSm1 = zeros(Ny,Nx,K);
uGSm1(1,:,:) = 0;
gm1 = 0;
uGSm1(Ny,:,:) = 100;
uGSm1(:,Nx,:) = 100;

uGSm1Stored = zeros(Ny,Nx,K);

for k = 1:K;
    for i = 2:Nx-1;
        for j = 2:Ny-1;
            if i == 2;
               uGSm1(j,i) = 1/3*(uGSm1(j-1,i) + uGSm1(j,i+1) + uGSm1(j+1,i) + deltax*gm1);
               else
               uGSm1(j,i) = 1/4*(uGSm1(j-1,i) + uGSm1(j,i+1) + uGSm1(j+1,i) + uGSm1(j,i-1));
            end
        end
    end
    uGSm1Stored(:,:,k) = uGSm1(:,:,1);
end

C = K;

for c = 1:C-1
    GSL2m1(c) = sum(sqrt(sum((uGSm1Stored(:,:,c+1) - uGSm1Stored(:,:,c))^2)));
end

KK = 101;

uGSm2 = zeros(Ny,Nx,KK);
uGSm2(1,:,:) = 0;
gm2 = 0;
uGSm2(Ny,:,:) = 100;
uGSm2(:,Nx,:) = 100;

uGSm2Stored = zeros(Ny,Nx,KK);

for kk = 1:KK;
    for ii = 2:Nx-1;
        for jj = 2:Ny-1;
            if ii == 2;
               uGSm2(jj,ii) = 1/4*(uGSm2(jj-1,ii) + uGSm2(jj+1,ii) + 2*uGSm2(jj,ii+1) + 2*gm2);
            else
               uGSm2(jj,ii) = 1/4*(uGSm2(jj-1,ii) + uGSm2(jj,ii+1) + uGSm2(jj+1,ii) + uGSm2(jj,ii-1));
            end
        end
    end
    uGSm2Stored(:,:,kk) = uGSm2(:,:,1);
end

CC = KK;

for cc = 1:CC-1
    GSL2m2(cc) = sum(sqrt(sum((uGSm2Stored(:,:,cc+1) - uGSm2Stored(:,:,cc))^2)));
end

%%%%% PLOTTING %%%%%

% plot(GSL2m1,'*r');
% hold on;
% plot (GSL2m2,'*b');
% title('L2 Norm of Residual');
% legend('1st order','2nd order');
% ylabel('L2 Norm');
% xlabel('Number of Iterations');
% hold off;
% 
% pause

% mesh(uGSm1(:,:,1));
% title('10x10 Solution with 1st order Method');
% ylabel('Ny');
% xlabel('Nx');
% zlabel('u');
% pause

% mesh(uGSm2(:,:,1));
% title('10x10 Solution with 2nd order Method');
% ylabel('Ny');
% xlabel('Nx');
% zlabel('u');

y = 0:deltay:1;

% plot(y,uGSm1(:,16,1),'*r');
% hold on;
% plot(Gauss20x20neumann(:,1),Gauss20x20neumann(:,2),'*b');
% set(gca,'FontSize',18);
% title('Solution of T with approx. 900 cells','FontSize',18);
% legend('MATLAB (finite difference)','OpenFOAM (finite volume)');
% xlabel('y','FontSize',18);
% ylabel('T(y)','FontSize',18);

plot(y,uGSm1(:,31,1),'*r');
hold on;
plot(Gauss40x40neumann(:,1),Gauss40x40neumann(:,2),'*b');
set(gca,'FontSize',18);
title('Solution of T with approx. 3500 cells','FontSize',18);
legend('MATLAB (finite difference)','OpenFOAM (finite volume)');
xlabel('y','FontSize',18);
ylabel('T(y)','FontSize',18);

% pause

% plot(Gauss40x40neumann(:,1),Gauss40x40neumann(:,3),'*r');
% hold on;
% plot(LeastSquares40x40neumann(:,1),LeastSquares40x40neumann(:,3),'*b');
% set(gca,'FontSize',18);
% title('Solution of grad(T) with approx. 3500 cells','FontSize',18);
% legend('Gauss','Least Squares');
% xlabel('y','FontSize',18);
% ylabel('grad(T(y))','FontSize',18);
% hold off;
% 
% pause
% 
% plot(Gauss20x20neumann(:,1),Gauss20x20neumann(:,3),'*r');
% hold on;
% plot(LeastSquares20x20neumann(:,1),LeastSquares20x20neumann(:,3),'*b');
% set(gca,'FontSize',18);
% title('Solution of grad(T) with approx. 900 cells','FontSize',18);
% legend('Gauss','Least Squares');
% xlabel('y','FontSize',18);
% ylabel('grad(T(y))','FontSize',18);
% hold off;