clear all;
format long;

load Gauss20x20dirichlet;
load LeastSquares20x20dirichlet;
load Gauss40x40dirichlet;
load LeastSquares40x40dirichlet;

a = 1;
b = 1;

Nx = 61;
Ny = 61;

deltax = a/(Nx-1);
deltay = b/(Ny-1);

K = 1000;

uJacobi = zeros(Ny,Nx,K+1);
uJacobi(1,:,:) = 0;
uJacobi(:,1,:) = 0;
uJacobi(Nx,:,:) = 0;
uJacobi(:,Ny,:) = 100;

for k = 1:K;
    for i = 2:Nx-1;
        for j = 2:Ny-1;
            uJacobi(j,i,k+1) = 1/4*(uJacobi(j-1,i,k) + uJacobi(j,i+1,k) + uJacobi(j+1,i,k) + uJacobi(j,i-1,k));
        end
    end
    JacobiL2(k) = sum(sqrt(sum((uJacobi(:,:,k+1) - uJacobi(:,:,k))^2)));
end

% surf(uJacobi(:,:,K+1));
% pause;

KK = 5000;

uGS = zeros(Ny,Nx,KK);
uGS(1,:,:) = 0;
uGS(:,1,:) = 0;
uGS(Nx,:,:) = 0;
uGS(:,Ny,:) = 100;

uGSStored = zeros(Ny,Nx,KK);

for kk = 1:KK;
    for ii = 2:Nx-1;
        for jj = 2:Ny-1;
            uGS(jj,ii) = 1/4*(uGS(jj-1,ii) + uGS(jj,ii+1) + uGS(jj+1,ii) + uGS(jj,ii-1));
        end
    end
    uGSStored(:,:,kk) = uGS(:,:,1);
end

CC = KK;

for cc = 1:CC-1
    GSL2(cc) = sum(sqrt(sum((uGSStored(:,:,cc+1) - uGSStored(:,:,cc))^2)));
end

% surf(uGS(:,:,1));

% uAnalytical = zeros(Ny,Nx);
% 
% for n = 1:2:201;
%     for x = 1:Nx-2;
%         for y = 1:Ny-2;
%             uAnalytical(y,x) = uAnalytical(y,x) + (1/n)*(sinh((n*pi*x*deltax)/b)/sinh((n*pi*a)/b))*sin((n*pi*y*deltay)/b);
%         end
%     end
% end
% 
% uAnalytical = (400/pi)*uAnalytical;
% uAnalytical = circshift(uAnalytical,[1 1]);

%%%%% PLOTTING %%%%%

% plot(GSL2,'*r');
% hold on;
% plot (JacobiL2,'*b');
% title('L2 Norm of Residual');
% legend('Gauss-Seidel','Jacobi');
% ylabel('L2 Norm');
% xlabel('Number of Iterations');
% hold off;
% 
% pause
% 
% mesh(uGS(:,:,1))
% 
% pause;

y = 0:deltay:1;

% plot(y,uGS(16,:,1),'*r');
% hold on;
% plot(Gauss20x20dirichlet(:,1),Gauss20x20dirichlet(:,2),'*b');
% set(gca,'FontSize',18);
% title('Solution of with approx. 900 cells','FontSize',18);
% legend('MATLAB (finite difference)','OpenFOAM (finite volume)');
% xlabel('y','FontSize',18);
% ylabel('T(y)','FontSize',18);

plot(y,uGS(31,:,1),'*r');
hold on;
plot(Gauss40x40dirichlet(:,1),Gauss40x40dirichlet(:,2),'*b');
set(gca,'FontSize',18);
title('Solution of T with approx. 3500 cells','FontSize',18);
legend('MATLAB (finite difference)','OpenFOAM (finite volume)');
xlabel('y','FontSize',18);
ylabel('T(y)','FontSize',18);

% pause

% plot(Gauss40x40dirichlet(:,1),Gauss40x40dirichlet(:,4),'*r');
% hold on;
% plot(LeastSquares40x40dirichlet(:,1),LeastSquares40x40dirichlet(:,4),'*b');
% set(gca,'FontSize',18);
% title('Solution of grad(T) with approx. 3500 cells','FontSize',18);
% legend('Gauss','Least Squares');
% xlabel('y','FontSize',18);
% ylabel('grad(T(y))','FontSize',18);
% hold off;
% 
% pause
% 
% plot(Gauss20x20dirichlet(:,1),Gauss20x20dirichlet(:,4),'*r');
% hold on;
% plot(LeastSquares20x20dirichlet(:,1),LeastSquares20x20dirichlet(:,4),'*b');
% set(gca,'FontSize',18);
% title('Solution of grad(T) with approx. 900 cells','FontSize',18);
% legend('Gauss','Least Squares');
% xlabel('y','FontSize',18);
% ylabel('grad(T(y))','FontSize',18);
% hold off;