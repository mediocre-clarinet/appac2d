function rightRunningWakes = solveWake(foils,Ainv,RHS,CT,options)
% Default options
opts.MaxIterations = 50;
opts.FunctionTolerance = 1e-6;
opts.RelaxationFactor = 0.5;
opts.NumPanels = 201;
opts.WakeLengthChords = 7;
opts.GammaDecayChords = 2;
opts.Display = 'iter';

% Attach flat wake to the trailing edge of each propulsive element
rightRunningWakes.m = opts.NumPanels + [0 0];
for i = 2:-1:1
    k = (i-1)*opts.NumPanels+(1:opts.NumPanels);
    rightRunningWakes.xo(k,:) = foils.xo((i-1)*foils.m(1)+1) + linspace(0,opts.WakeLengthChords,opts.NumPanels).';
    rightRunningWakes.yo(k,:) = foils.yo((i-1)*foils.m(1)+1) + zeros(opts.NumPanels,1);
    rightRunningWakes.dx(k,:) = [diff(rightRunningWakes.xo(k)); 1e3];
    rightRunningWakes.dy(k,:) = zeros(opts.NumPanels,1);
end
rightRunningWakes.theta = atan2(rightRunningWakes.dy,rightRunningWakes.dx);
collocPoints = [rightRunningWakes.xo+rightRunningWakes.dx/2 rightRunningWakes.yo+rightRunningWakes.dy/2];
leftRunningWakes.m = rightRunningWakes.m;
leftRunningWakes.xo = rightRunningWakes.xo + rightRunningWakes.dx;
leftRunningWakes.dx = -rightRunningWakes.dx;

gammaInf = sqrt(2*CT + 1) - 1;
gamma = repelem([gammaInf;-gammaInf],opts.NumPanels+1);

for iter = 1:30
    leftRunningWakes.yo = rightRunningWakes.yo + rightRunningWakes.dy;
    leftRunningWakes.dy = -rightRunningWakes.dy;
    leftRunningWakes.theta = atan2(leftRunningWakes.dy,leftRunningWakes.dx);

    [u,v] = influence(foils.co,rightRunningWakes);
    A = -u.*sin(foils.theta) + v.*cos(foils.theta);
    g = Ainv*(RHS - [A*gamma;-gamma(1);-gamma(opts.NumPanels+1);zeros(numel(foils.m)-2,1)]);

    [uu,vv] = influence(collocPoints,foils);
    u = uu*g + 1; v = vv*g;
    [uu,vv] = influence(collocPoints,rightRunningWakes);
    u = u + uu/2*gamma; v = v + vv/2*gamma;
    [uu,vv] = influence(collocPoints,leftRunningWakes);
    u = u + uu/2*gamma; v = v + vv/2*gamma;

    rightRunningWakes.dy = v./u.*rightRunningWakes.dx;
    rightRunningWakes.dy(opts.NumPanels) = 0;
    rightRunningWakes.dy(end) = 0;
    rightRunningWakes.yo(2:opts.NumPanels) = rightRunningWakes.yo(1) + cumsum(rightRunningWakes.dy(1:opts.NumPanels-1));
    rightRunningWakes.yo(opts.NumPanels+2:end) = rightRunningWakes.yo(opts.NumPanels+1) + cumsum(rightRunningWakes.dy(opts.NumPanels+1:end-1));
    rightRunningWakes.theta = atan2(rightRunningWakes.dy,rightRunningWakes.dx);
    collocPoints(:,2) = rightRunningWakes.yo + rightRunningWakes.dy/2;

    gamma([1:opts.NumPanels opts.NumPanels+2:end-1]) = CT./(u.^2 + v.^2);
    gamma(opts.NumPanels+2:end-1) = -gamma(opts.NumPanels+2:end-1);

%    figure
%    plot(0:opts.NumPanels,reshape(gamma,[],2))
end
