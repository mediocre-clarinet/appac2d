function [u,v] = influence(co,af,part)
% INFLUENCE  Calculate aerodynamic influence of linear-strength vortex panels.
%   [U,V] = INFLUENCE(CO,AF,PART) returns coefficient matrices U and V, which
%   when multiplied by the circulation values associated with the surfaces in
%   the unified struct AF, yield the x- and y-components of velocity induced
%   at the collocation points CO, respectively. CO is an N-by-2 array of 2D
%   coordinates. If a collocation point lies on a panel, PART determines the
%   side from which the solution is taken (top:1, bottom:-1, averaged:0).

% Convert control point to local panel coords
xt = co(:,1) - af.xo.';
yt = co(:,2) - af.yo.';

% Precompute trig expressions as 1-by-m arrays
costh = cos(af.theta).';
sinth = sin(af.theta).';

% Find control point coords in CS of all panels
xp = xt.*costh + yt.*sinth;
yp = -xt.*sinth + yt.*costh;
x2 = af.dx.'.*costh + af.dy.'.*sinth; % 1-by-m

% Find theta1, theta2, r1, r2
theta1 = atan2(yp,xp);
theta2 = atan2(yp,xp-x2);
dtheta = theta2 - theta1; % precompute
dtheta((abs(yp) < 1e-12) & (xp > 0) & (xp < x2)) = part*pi;

%r1 = sqrt(xp.^2 + yp.^2); r2 = sqrt((xp-x2).^2 + yp.^2); ln = log(r2./r1);
ln = 0.5*log(((xp-x2).^2 + yp.^2)./(xp.^2 + yp.^2)); % precompute

% Compute influence coefficients
ap = yp.*ln + xp.*dtheta;
bp = xp.*ln + x2 - yp.*dtheta;
am = -yp.*ln + (x2 - xp).*dtheta;
bm = (x2 - xp).*ln - x2 + yp.*dtheta;
c = 1./(2*pi*x2);

% Velocity decomposition used in Katz and Plotkin Eq. 11.103
ua = c.*(am.*costh - bm.*sinth);
ub = c.*(ap.*costh - bp.*sinth);
va = c.*(am.*sinth + bm.*costh);
vb = c.*(ap.*sinth + bp.*costh);

i = cumsum([af.m]+1); % columns where zeros are written
j = 1:i(end); j(i) = []; % a-shifted column indices
u(:,j+1) = ub; % write to b-shifted columns
u(:,j) = u(:,j) + ua;
v(:,j+1) = vb;
v(:,j) = v(:,j) + va;