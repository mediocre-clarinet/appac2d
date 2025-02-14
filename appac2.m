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
Ainv = inv(A);
RHS = [sin(theta);zeros(nSurfs,1)];
gammaInf = sqrt(2*CT + 1) - 1;

wakes(1).gamma = zeros(N+1,1) + gammaInf;
wakes(2).gamma = zeros(N+1,1) - gammaInf;

AA = zeros(M,2*(N+1));
BB = zeros(M,2*(N+1));

for iter = 1:3
    for i = 1:2
        [u,v] = influence(co,wakes(i));
        AA(:,(i-1)*(N+1)+(1:N+1)) = -u.*sin(theta) + v.*cos(theta);
        BB(:,(i-1)*(N+1)+(1:N+1)) = u.*cos(theta) + v.*sin(theta);
    end
    
    augmentorRHS = [AA*vertcat(wakes.gamma);-wakes(1).gamma(1);-wakes(2).gamma(end)];
    
    gamma = Ainv*(RHS - augmentorRHS);
    tmp = mat2cell(gamma,[foils.m]+1);
    [foils(:).gamma] = tmp{:};
    
    wco = vertcat(wakes.co);
    u = 1; v = 0;
    for i = 1:2
        [uu,vv] = influence(wco,wakes(i));
        u = u + uu*wakes(i).gamma;
        v = v + vv*wakes(i).gamma;
        [uu,vv] = influence(wco,foils(i));
        u = u + uu*foils(i).gamma;
        v = v + vv*foils(i).gamma;
    end
    
    wakes(1).dy = v(1:N)./u(1:N).*wakes(1).dx;
    wakes(1).dy(N) = 0;
    wakes(2).dy = v(N+1:2*N)./u(N+1:2*N).*wakes(2).dx;
    wakes(2).dy(1) = 0;
    wakes(1).y(2:N+1) = wakes(1).y(1) + cumsum(wakes(1).dy);
    wakes(2).y(2:N+1) = wakes(2).y(1) + cumsum(wakes(2).dy);
    wakes(2).y = wakes(2).y - wakes(2).y(N+1) + foils(2).y(1);
    wakes(1).theta = atan2(wakes(1).dy,wakes(1).dx);
    wakes(2).theta = atan2(wakes(2).dy,wakes(2).dx);
    wakes(1).gamma(1:N) = CT./(u(1:N).^2 + v(1:N).^2);
    wakes(1).gamma(N:N+1) = gammaInf;
    wakes(2).gamma(2:N+1) = CT./(u(N+1:2*N).^2 + v(N+1:2*N).^2);
    wakes(2).gamma(1:2) = -gammaInf;
    wakes(1).gamma(41:N) = linspace(wakes(1).gamma(41),wakes(1).gamma(N),161);
    wakes(2).gamma(2:162) = linspace(wakes(2).gamma(2),wakes(2).gamma(162),161);

    figure
    plot(u(1:N))
    hold on
    plot(flipud(u(N+1:2*N)))
end

Qtan = B*gamma + cos(theta);
Cp = 1 - Qtan.^2;
xc = co*R(1,:).';