clear all;
format long;

load stepU12fine;
load stepU21fine;
load stepU30fine;
load stepU39fine;
load stepU48fine;
load stepU57fine;
load stepU66fine;
load stepU75fine;
load stepU84fine;
load stepU93fine;
load stepU102fine;
load stepU111fine;
load stepU120fine;
load stepU129fine;
load stepU138fine;
load stepU147fine;
load stepU156fine;
load stepU291fine;
load stepU300fine;

load x7U12fine;
load x7U21fine;
load x7U30fine;
load x7U39fine;
load x7U48fine;
load x7U57fine;
load x7U66fine;
load x7U75fine;
load x7U84fine;
load x7U93fine;
load x7U102fine;
load x7U111fine;
load x7U120fine;
load x7U129fine;
load x7U138fine;
load x7U147fine;
load x7U156fine;
load x7U291fine;
load x7U300fine;

load stepU300medium;
load x7U300medium;

load stepU300coarse;
load x7U300coarse;
load x7U300coarseLinear;
load x7U300coarseUpwind;

%%%%% GARTLING'S EXPERIMENT %%%%%

y = 0:0.001:0.5;

stepUGartling = 24*y.*(0.5 - y);

%%%%% PLOTTING %%%%%

%%%%% No. 1 %%%%%

plot(x7U12fine(:,2),x7U12fine(:,1),'ro-');
set(gca,'FontSize',18);
hold on;
plot(x7U21fine(:,2),x7U21fine(:,1),'bo-');
plot(x7U30fine(:,2),x7U30fine(:,1),'go-');
plot(x7U39fine(:,2),x7U39fine(:,1),'ko-');
plot(x7U48fine(:,2),x7U48fine(:,1),'mo-');
plot(x7U57fine(:,2),x7U57fine(:,1),'yo-');
plot(x7U66fine(:,2),x7U66fine(:,1),'r*-');
plot(x7U75fine(:,2),x7U75fine(:,1),'b*-');
% plot(x7U84fine(:,2),x7U84fine(:,1),'g*-');
% plot(x7U93fine(:,2),x7U93fine(:,1),'k*-');
% plot(x7U102fine(:,2),x7U102fine(:,1),'m*-');
% plot(x7U111fine(:,2),x7U111fine(:,1),'y*-');
title('Horizontal Velocity at x=7 for Various Times','FontSize',18);
legend('t=12 sec','t=21 sec','t=30 sec','t=39 sec','t=48 sec','t=57 sec','t=66 sec','t=75 sec');
xlabel('horizontal velocity','FontSize',18);
ylabel('y','FontSize',18);
hold off;
pause;

plot(x7U12fine(:,3),x7U12fine(:,1),'ro-');
set(gca,'FontSize',18);
hold on;
plot(x7U21fine(:,3),x7U21fine(:,1),'bo-');
plot(x7U30fine(:,3),x7U30fine(:,1),'go-');
plot(x7U39fine(:,3),x7U39fine(:,1),'ko-');
plot(x7U48fine(:,3),x7U48fine(:,1),'mo-');
plot(x7U57fine(:,3),x7U57fine(:,1),'yo-');
plot(x7U66fine(:,3),x7U66fine(:,1),'r*-');
plot(x7U75fine(:,3),x7U75fine(:,1),'b*-');
plot(x7U84fine(:,3),x7U84fine(:,1),'g*-');
plot(x7U93fine(:,3),x7U93fine(:,1),'k*-');
plot(x7U102fine(:,3),x7U102fine(:,1),'m*-');
plot(x7U111fine(:,3),x7U111fine(:,1),'y*-');
plot(x7U120fine(:,3),x7U120fine(:,1),'r^-');
plot(x7U129fine(:,3),x7U129fine(:,1),'b^-');
plot(x7U138fine(:,3),x7U138fine(:,1),'g^-');
plot(x7U147fine(:,3),x7U147fine(:,1),'k^-');
plot(x7U291fine(:,3),x7U291fine(:,1),'m^-');
plot(x7U300fine(:,3),x7U300fine(:,1),'y^-');
title('Vertical Velocity at x=7 for Various Times','FontSize',18);
legend('t=12 sec','t=21 sec','t=30 sec','t=39 sec','t=48 sec','t=57 sec','t=66 sec','t=75 sec','t=84 sec','t=93 sec','t=102 sec','t=111 sec','t=120 sec','t=129 sec','t=138 sec','t=147 sec','t=291 sec','t=300 sec');
xlabel('vertical velocity','FontSize',18);
ylabel('y','FontSize',18);
hold off;
pause;

%%%%% No. 2 %%%%%

plot(x7U300fine(:,2),x7U300fine(:,1),'ro-');
set(gca,'FontSize',18);
hold on;
plot(x7U300medium(:,2),x7U300medium(:,1),'bo-');
plot(x7U300coarse(:,2),x7U300coarse(:,1),'go-');
title('Horizontal Velocity at x=7 for Various Grids','FontSize',18);
xlabel('horizontal velocity','FontSize',18);
ylabel('y','FontSize',18);
legend('fine (54,999 cells)','medium (27,417 cells)','coarse (13,221 cells)');
hold off;
pause;

%%%%% No. 3 %%%%%

plot(x7U300coarse(:,2),x7U300coarse(:,1),'ro-');
set(gca,'FontSize',18);
hold on;
plot(x7U300coarseUpwind(:,2),x7U300coarseUpwind(:,1),'bo-');
plot(x7U300coarseLinear(:,2),x7U300coarseLinear(:,1),'go-');
title('Horizontal Velocity at x=7 for Various Divergence Schemes','FontSize',18);
xlabel('horizontal velocity','FontSize',18);
ylabel('y','FontSize',18);
legend('Linear Upwind','Upwind','Linear');
hold off;
pause;

%%%%% No. 4 %%%%%

plot(x7U300fine(:,2),x7U300fine(:,1),'ro');
set(gca,'FontSize',18);
title('Horizontal Velocity at x=7','FontSize',18);
xlabel('horizontal velocity','FontSize',18);
ylabel('y','FontSize',18);
pause;

plot(x7U300fine(:,3),x7U300fine(:,1),'ro');
set(gca,'FontSize',18);
title('Vertical Velocity at x=7','FontSize',18);
xlabel('vertical velocity','FontSize',18);
ylabel('y','FontSize',18);
pause;

plot(stepU300fine(:,2),stepU300fine(:,1),'ro');
hold on;
plot(stepUGartling,y,'b-');
title('Horizontal Velocity at Step','FontSize',18);
xlabel('horizontal velocity','FontSize',18);
ylabel('y','FontSize',18);
legend('OpenFOAM','Gartling');
hold off;
pause;
