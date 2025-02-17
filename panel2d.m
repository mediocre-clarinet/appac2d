function [foils,wakes,Cp,xc] = panel2d(surfaces,alphaDeg,xDisk,CT)
% PANEL2D  Panel method in two dimensions.
nSurfs = numel(surfaces);
R = [cosd(alphaDeg) -sind(alphaDeg);sind(alphaDeg) cosd(alphaDeg)];

% Determine total number of surface panels for memory allocation %%%%%%%%%%%%%
for i = nSurfs:-1:1 % populating the last entry allocates necessary space
    foils.m(i) = size(surfaces{i},1) - 1;
end
M = sum([foils.m]); % total number of panels

% Build unified struct for all lifting surfaces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = M;
for i = nSurfs:-1:1
    k = k - foils.m(i); % shift indices to rows in array for writing
    coords = surfaces{i}*R; % rotate coordinates by AoA
    foils.xo(1+k:foils.m(i)+k,:) = coords(1:end-1,1);
    foils.yo(1+k:foils.m(i)+k,:) = coords(1:end-1,2);
    foils.dx(1+k:foils.m(i)+k,:) = diff(coords(:,1));
    foils.dy(1+k:foils.m(i)+k,:) = diff(coords(:,2));
end
foils.theta = atan2(foils.dy,foils.dx);
foils.co = [foils.xo+foils.dx/2 foils.yo+foils.dy/2];

% Create aerodynamic influence coefficient matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(M+nSurfs);
for i = 1:nSurfs % Kutta condition
    A(M+i,k+1) = 1;
    A(M+i,k+foils.m(i)+1) = 1;
    k = k + foils.m(i) + 1;
end
[U,V] = influence(foils.co,foils,1);
A(1:M,:) = -U.*sin(foils.theta) + V.*cos(foils.theta); % normal component
B = U.*cos(foils.theta) + V.*sin(foils.theta); % tangent component
RHS = [sin(foils.theta);zeros(nSurfs,1)];

% Solve for global circulation solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wakes,foils.gamma,iter,E] = solveWake(foils,inv(A),RHS,CT);

[U,V] = influence(foils.co,wakes,1);
D = U.*cos(foils.theta) + V.*sin(foils.theta);

Qtan = B*foils.gamma + cos(foils.theta) + D*wakes.gamma;
Cp = 1 - Qtan.^2;
xc = foils.co*R(1,:).';

% Correct Cp aft of the actuator disk where the total pressure is higher %%%%%
k1 = find(xc(1:foils.m(1)) < xDisk, 1, 'last');
k2 = find(xc(foils.m(1)+(1:foils.m(2))) < xDisk, 1, 'first');
Cp(k1+1:foils.m(1)+k2-1) = Cp(k1+1:foils.m(1)+k2-1) + 2*CT;

% Contour plotting under construction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = wakes.m(1);
% Distinguish jet and outer flow evaluation by making one wake left-running
wakes.xo(1:N) = wakes.xo(1:N) + wakes.dx(1:N);
wakes.yo(1:N) = wakes.yo(1:N) + wakes.dy(1:N);
wakes.dx(1:N) = -wakes.dx(1:N);
wakes.dy(1:N) = -wakes.dy(1:N);
wakes.theta(1:N) = atan2(wakes.dy(1:N),wakes.dx(1:N));

node = [foils.co(k1+1:foils.m(1),:); ...
        wakes.co(1:N-1,:); ...
        flipud(wakes.co(N+1:2*N-1,:)); ...
        foils.co(foils.m(1)+(1:k2-1),:)];

edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;
opts.kind = 'delfront';
opts.rho2 = 1;
[vert,etri,tria,tnum] = refine2(node,edge,[],opts);

[U,V] = influence(vert,foils,1);
u = U*foils.gamma + 1;
v = V*foils.gamma;
[U,V] = influence(vert,wakes,-1); % -1 for jet interior
u = u + U*wakes.gamma;
v = v + V*wakes.gamma;

figure;
cmap = crameri('berlin');
colormap(cmap);
patch('Faces',tria(:,1:3),'Vertices',vert, ...
      'FaceVertexCData',sqrt(u.*u + v.*v), ...
      'FaceColor','interp','EdgeColor','none');
hold on; axis image off;

node = [flipud(wakes.co(1:N-1,:)); ...
        foils.co(1:k1+1,:); ...
        foils.co(foils.m(1)+(k2-1:foils.m(2)),:); ...
        wakes.co(N+1:2*N-1,:); ...
        wakes.co(2*N-1,1) 2; ...
        -2 2;-2 -2; ...
        wakes.co(N-1,1) -2];

edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;
opts.kind = 'delfront';
opts.rho2 = 1;
[vert,etri,tria,tnum] = refine2(node,edge,[],opts);

[U,V] = influence(vert,foils,1);
u = U*foils.gamma + 1;
v = V*foils.gamma;
[U,V] = influence(vert,wakes,1);
u = u + U*wakes.gamma;
v = v + V*wakes.gamma;

patch('Faces',tria(:,1:3),'Vertices',vert, ...
      'FaceVertexCData',sqrt(u.*u + v.*v), ...
      'FaceColor','interp','EdgeColor','none');

cl = get(gca,'CLim');
set(gca,'CLim',max(abs(cl-1))*[-1 1]+1);