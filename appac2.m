function [foils,surfaceVortexSheets,wakes,Cp,xc] = appac2(surfaces,alphaDeg,xDisk,CT,options)
nSurfs = numel(surfaces);
R = [cosd(alphaDeg) -sind(alphaDeg);sind(alphaDeg) cosd(alphaDeg)];
k = find(surfaces{1}(:,1) <= xDisk, 1, 'last' );
surfaceVortexSheets(1) = initSurface([surfaces{1}(k,:) + (surfaces{1}(k+1,:) - surfaces{1}(k,:))*(xDisk - surfaces{1}(k,1))/(surfaces{1}(k+1,1) - surfaces{1}(k,1)); surfaces{1}(k+1:end,:)]*R);
k = find(surfaces{2}(:,1) <= xDisk, 1, 'first');
surfaceVortexSheets(2) = initSurface([surfaces{2}(1:k-1,:); surfaces{2}(k-1,:) + (surfaces{2}(k,:) - surfaces{2}(k-1,:))*(xDisk - surfaces{2}(k-1,1))/(surfaces{2}(k,1) - surfaces{2}(k-1,1))]*R);

% Determine the total number of surface panels so that we can preallocate space for variables
M = 0;
for i = nSurfs:-1:1 % populating an array from the last entry first allocates necessary space
    foils.m(i) = size(surfaces{i},1) - 1;
    M = M + foils.m(i);
end

k = M;
for i = nSurfs:-1:1
    k = k - foils.m(i); % shift indices to rows in array for writing
    coords = surfaces{i}*R; % rotate coordinates by AoA
    foils.xo(k+(1:foils.m(i)),:) = coords(1:end-1,1);
    foils.yo(k+(1:foils.m(i)),:) = coords(1:end-1,2);
    foils.dx(k+(1:foils.m(i)),:) = diff(coords(:,1));
    foils.dy(k+(1:foils.m(i)),:) = diff(coords(:,2));
end
foils.theta = atan2(foils.dy,foils.dx);
foils.co = [foils.xo+foils.dx/2 foils.yo+foils.dy/2];

% Create aerodynamic influence coefficient matrix
A = zeros(M+nSurfs);
[u,v] = influence(foils.co,foils,pi);
A(1:M,:) = -u.*sin(foils.theta) + v.*cos(foils.theta); % normal component
B = u.*cos(foils.theta) + v.*sin(foils.theta); % tangent component
% Kutta condition
k = 1;
for i = 1:nSurfs
    A(M+i,k) = 1;
    A(M+i,k+foils.m(i)) = 1;
    k = k + foils.m(i) + 1;
end
RHS = [sin(foils.theta);zeros(nSurfs,1)];

wakes = solveWake(foils,inv(A),RHS,CT,options);

gamma = A \ RHS;

Qtan = B*gamma + cos(foils.theta);
Cp = 1 - Qtan.^2;
xc = foils.co*R(1,:).';
end
