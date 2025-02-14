function [foils,surfaceVortexSheets,wakes,Cp,xc] = appac2(surfaces,alphaDeg,xDisk,CT,options)
nSurfs = numel(surfaces);
R = [cosd(alphaDeg) -sind(alphaDeg);sind(alphaDeg) cosd(alphaDeg)];
k = find(surfaces{1}(:,1) <= xDisk, 1, 'last' );
surfaceVortexSheets(1) = initSurface([surfaces{1}(k,:) + (surfaces{1}(k+1,:) - surfaces{1}(k,:))*(xDisk - surfaces{1}(k,1))/(surfaces{1}(k+1,1) - surfaces{1}(k,1)); surfaces{1}(k+1:end,:)]*R);
k = find(surfaces{2}(:,1) <= xDisk, 1, 'first');
surfaceVortexSheets(2) = initSurface([surfaces{2}(1:k-1,:); surfaces{2}(k-1,:) + (surfaces{2}(k,:) - surfaces{2}(k-1,:))*(xDisk - surfaces{2}(k-1,1))/(surfaces{2}(k,1) - surfaces{2}(k-1,1))]*R);
for i = nSurfs:-1:1
    foils(i) = initSurface(surfaces{i}*R);
end
N = 201;
wakes(1) = initSurface([foils(1).x(end)+[linspace(0,10,N) 1010].' foils(1).y(end)+zeros(N+1,1)]);
wakes(2) = initSurface([foils(2).x(1)+[1010 linspace(10,0,N)].' foils(2).y(1)+zeros(N+1,1)]);

M = sum([foils.m]);
A = zeros(M+nSurfs);
B = zeros(M,M+nSurfs);

co = vertcat(foils.co);
theta = vertcat(foils.theta);
k = 1;
for i = 1:nSurfs
    [u,v] = influence(co,foils(i));
    A(1:M,k+(0:foils(i).m)) = -u.*sin(theta) + v.*cos(theta);
    B(:,k+(0:foils(i).m)) = u.*cos(theta) + v.*sin(theta);
    A(M+i,k) = 1;
    A(M+i,k+foils(i).m) = 1;
    k = k + foils(i).m + 1;
end
RHS = [sin(theta);zeros(nSurfs,1)];
gamma = A \ RHS;
Qtan = B*gamma + cos(theta);
Cp = 1 - Qtan.^2;
xc = co*R(1,:).';