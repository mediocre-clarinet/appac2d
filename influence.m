function [u,v] = influence(co,af)
% co: N-by-2 array of collocation points where influences are calculated
% af: initialized surface
n = size(co,1);
m = af.m; % number of panels

% convert control point to local panel coords
xt = co(:,1) - af.x(1:m).';
yt = co(:,2) - af.y(1:m).';

% precompute trig expressions as 1-by-m arrays
costh = cos(af.theta).';
sinth = sin(af.theta).';

% find control point coords in CS of all panels
xp = xt.*costh + yt.*sinth;
yp = -xt.*sinth + yt.*costh;
x2 = af.dx.'.*costh + af.dy.'.*sinth; % 1-by-m
%y2 = 0;

% find theta1, theta2, r1, r2
theta1 = atan2(yp,xp);
theta2 = atan2(yp,xp-x2);
dtheta = theta2 - theta1; % precompute
k = (abs(yp) < 1e-12) & (xp > 0) & (xp < x2);
dtheta(k) = pi; % fix angular difference for self-induction

%r1 = sqrt(xp.^2 + yp.^2);
%r2 = sqrt((xp-x2).^2 + yp.^2);
%ln = log(r2./r1); % precompute
ln = 0.5*log(((xp-x2).^2 + yp.^2)./(xp.^2 + yp.^2)); % precompute

% compute influence coefficients
ap = yp.*ln + xp.*dtheta;
bp = xp.*ln + x2 - yp.*dtheta;
am = -yp.*ln + (x2 - xp).*dtheta;
bm = (x2 - xp).*ln - x2 + yp.*dtheta;
c = 1./(2*pi*x2);

ua = c.*(am.*costh - bm.*sinth);
ub = c.*(ap.*costh - bp.*sinth);
va = c.*(am.*sinth + bm.*costh);
vb = c.*(ap.*sinth + bp.*costh);

u = [ua, zeros(n,1)] + [zeros(n,1), ub];
v = [va, zeros(n,1)] + [zeros(n,1), vb];
end