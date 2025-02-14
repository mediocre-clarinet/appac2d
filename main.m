clear
close all

fid = fopen('mainVec0.dat','r');
coords = fscanf(fid,'%f %f',[2 Inf]);
fclose(fid);

af = initSurface(coords.');

[u,v] = influence(af.co,af);
A = -u.*sin(af.theta) + v.*cos(af.theta);
B = u.*cos(af.theta) + v.*sin(af.theta);
A(af.m+1,1) = 1; A(end) = 1;
RHS = [sin(af.theta);0];

gamma = A \ RHS;

Qtan = B*gamma + cos(af.theta);
Cp = 1 - Qtan.^2;
Cl = -Cp.'*af.dx;
Cd = Cp.'*af.dy;

figure
plot(af.co(:,1),Cp)
set(gca,'YDir','reverse')