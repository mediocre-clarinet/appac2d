function s = initSurface(coords)
% coords: N-by-2 array of connected coordinates
m = size(coords,1) - 1; % number of panels
s.m = m;

s.x = coords(:,1);
s.y = coords(:,2);

% find panel angles
s.dx = s.x(2:m+1) - s.x(1:m);
s.dy = s.y(2:m+1) - s.y(1:m);
s.theta = atan2(s.dy,s.dx);

% establish control points
%s.xc = 0.5*(s.x(1:m) + s.x(2:m+1));
%s.yc = 0.5*(s.y(1:m) + s.y(2:m+1));
s.co = 0.5*(coords(1:m,:) + coords(2:m+1,:));
end