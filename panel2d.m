function [Cp,xc] = panel2d(surfaces,alphaDeg,varargin)
% PANEL2D  Panel method in two dimensions.
%   PANEL2D(SURFACES,ALPHADEG)
%   PANEL2D(SURFACES,ALPHADEG,CT,XDISK)
%   PANEL2D(SURFACES,ALPHADEG,CT,XDISK,WAKEOPTIONS)
%   PANEL2D(___,NAME,VALUE)
[oper,CT,xDisk,wakeOptions,options] = parseInput(varargin);

nSurfs = numel(surfaces);

if ((oper == 2) && (nSurfs == 1)) || (nSurfs == 0)
    error('panel2d:invalidInput','Not enough surfaces.');
end

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
if oper == 1
    foils.gamma = A \ RHS;
    Qtan = B*foils.gamma + cos(foils.theta);
else
    [wakes,foils.gamma,iter,E] = solveWake(foils,inv(A),RHS,CT,wakeOptions);
    [U,V] = influence(foils.co,wakes,1);
    D = U.*cos(foils.theta) + V.*sin(foils.theta);
    Qtan = B*foils.gamma + cos(foils.theta) + D*wakes.gamma;
end

Cp = mat2cell(1 - Qtan.^2, foils.m);
xc = mat2cell(foils.co*R(1,:).', foils.m);

if oper == 2
    % Correct Cp aft of the actuator disk where the total pressure is higher
    k1 = find(xc{1} < xDisk, 1, 'last');
    k2 = find(xc{2} < xDisk, 1, 'first');
    Cp{1}(k1+1:end) = Cp{1}(k1+1:end) + 2*CT;
    Cp{2}(1:k2-1) = Cp{2}(1:k2-1) + 2*CT;
end

% Create contour plot of velocity magnitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meshing settings
opts.kind = 'delfront';
opts.rho2 = 1;

figure;
cmap = crameri(options.Colormap);
colormap(cmap);
hold on; axis image off;

if oper == 1
    j = cumsum(foils.m);
    node = [foils.co;4 2;-2 2;-2 -2;4 -2];
    edge = (1:j(end)+4).' + [0 1]; edge(j,2) = edge(j,2) - foils.m.';
    edge(end,2) = j(end) + 1;
    [vert,etri,tria,tnum] = refine2(node,edge,[],opts);
    [U,V] = influence(vert,foils,1);
    u = U*foils.gamma + 1;
    v = V*foils.gamma;
    patch('Faces',tria,'Vertices',vert, ...
        'FaceVertexCData',sqrt(u.*u + v.*v), ...
        'FaceColor','interp','EdgeColor','none');
    cl = get(gca,'CLim');
    set(gca,'CLim',max(abs(cl-1))*[-1 1]+1);
    FlowP = tristream(tria,vert(:,1),vert(:,2),u,v,zeros(1,41)-1.99,-1:0.05:1);
    for i = 1:numel(FlowP)
        plot(FlowP(i).x,FlowP(i).y,'-','Color',[.5 .5 .5]);
    end
    return
end

N = wakes.m(1);
% Distinguish jet and outer flow evaluation by making one wake left-running
wakes.xo(1:N) = wakes.xo(1:N) + wakes.dx(1:N);
wakes.yo(1:N) = wakes.yo(1:N) + wakes.dy(1:N);
wakes.dx(1:N) = -wakes.dx(1:N);
wakes.dy(1:N) = -wakes.dy(1:N);
wakes.theta(1:N) = atan2(wakes.dy(1:N),wakes.dx(1:N));

% Mesh jet volume
node = [foils.co(k1+1:foils.m(1),:); ...
        wakes.co(1:N-1,:); ...
        flipud(wakes.co(N+1:2*N-1,:)); ...
        foils.co(foils.m(1)+(1:k2-1),:)];

edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;

[vert{1},etri,tria{1},tnum] = refine2(node,edge,[],opts);

[U,V] = influence(vert{1},foils,1);
u{1} = U*foils.gamma + 1;
v{1} = V*foils.gamma;
[U,V] = influence(vert{1},wakes,-1); % -1 for jet interior
u{1} = u{1} + U*wakes.gamma;
v{1} = v{1} + V*wakes.gamma;

patch('Faces',tria{1},'Vertices',vert{1}, ...
    'FaceVertexCData',sqrt(u{1}.^2 + v{1}.^2), ...
    'FaceColor','interp','EdgeColor','none');

% Mesh outer volume
node = [flipud(wakes.co(1:N-1,:)); ...
        foils.co(1:k1+1,:); ...
        foils.co(foils.m(1)+(k2-1:foils.m(2)),:); ...
        wakes.co(N+1:2*N-1,:); ...
        wakes.co(2*N-1,1) 2; ...
        -2 2;-2 -2; ...
        wakes.co(N-1,1) -2];

edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;

if nSurfs > 2
    node2 = foils.co(foils.m(1)+foils.m(2)+1:M,:);
    edge2 = (1:size(node2,1)).' + [0 1] + size(node,1);
    j = cumsum(foils.m(3:nSurfs));
    edge2(j,2) = edge2(j,2) - foils.m(3:nSurfs).';
    node = [node;node2];
    edge = [edge;edge2];
end

[vert{2},etri,tria{2},tnum] = refine2(node,edge,[],opts);

[U,V] = influence(vert{2},foils,1);
u{2} = U*foils.gamma + 1;
v{2} = V*foils.gamma;
[U,V] = influence(vert{2},wakes,1);
u{2} = u{2} + U*wakes.gamma;
v{2} = v{2} + V*wakes.gamma;

patch('Faces',tria{2},'Vertices',vert{2}, ...
    'FaceVertexCData',sqrt(u{2}.^2 + v{2}.^2), ...
    'FaceColor','interp','EdgeColor','none');

cl = get(gca,'CLim');
set(gca,'CLim',max(abs(cl-1))*[-1 1]+1);

FlowP = tristream([tria{1};tria{2}+size(vert{1},1)], ...
            [vert{1}(:,1);vert{2}(:,1)],[vert{1}(:,2);vert{2}(:,2)], ...
            [u{1};u{2}],[v{1};v{2}],zeros(1,41)-1.99,-1:0.05:1);
for i = 1:numel(FlowP)
    plot(FlowP(i).x,FlowP(i).y,'-','Color',[.5 .5 .5]);
end
end


% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [oper,CT,xDisk,wakeOptions,options] = parseInput(args)
% PARSEINPUT  Handle input arguments from the various input syntaxes.
oper = 1; % operating mode (vanilla:1, appac:2)
CT = 0;   % values of CT and xDisk are inconsequential when oper=1
xDisk = 0;

% Default options
wakeOptions.MaxIterations = 50;
wakeOptions.FunctionTolerance = 1e-6;
wakeOptions.RelaxationFactor = 0.5;
wakeOptions.NumPanels = 100;
wakeOptions.WakeLengthChords = 9;
wakeOptions.GammaDecayChords = 3;
wakeOptions.NodeSpacing = 'cosine';
wakeOptions.Display = 'iter';

options.Colormap = 'berlin';
options.Mesh = 'off';
options.Plot = 'off';
options.CData = 'q';

% Perform checks on input arguments and parse accordingly %%%%%%%%%%%%%%%%%%%%
nArgs = numel(args);

if nArgs == 0
    return
end

if nArgs == 1
    error('panel2d:invalidInputSyntax','Invalid input syntax.');
end

if mod(nArgs,2) == 0
    if ~(isnumeric(args{1}) && isscalar(args{1}) && ...
         isnumeric(args{2}) && isscalar(args{2}))&& ...
        ~ischar   (args{1})
        error('panel2d:incorrectInputClass','Incorrect input class.');
    end
    NameValueStartIdx = isnumeric(args{1})*2 + 1;
else
    if (~isnumeric(args{1}) || ~isscalar(args{1}) || ...
        ~isnumeric(args{2}) || ~isscalar(args{2}) || ...
        ~isstruct (args{3}) )
        error('panel2d:incorrectInputClass','Incorrect input class.');
    end
    NameValueStartIdx = 4;
    wakeOptions = updateOptions(wakeOptions, ...
                    [fieldnames(args{3}) struct2cell(args{3})].');
end
if NameValueStartIdx > 2
    oper = 2;
    CT = args{1};
    xDisk = args{2};
end
options = updateOptions(options,args(NameValueStartIdx:nArgs));
end

function options = updateOptions(options,NameValuePairs)
% UPDATEOPTIONS  Update options using pairs of Name-Value arguments.
names = fieldnames(options);
for i = 1:2:numel(NameValuePairs)-1
    k = strcmpi(names,NameValuePairs{i});
    if ~any(k)
        error('panel2d:invalidOptions', ...
            sprintf('%s is not a recognized option.',NameValuePairs{i}));
    end
    options.(names{k}) = NameValuePairs{i+1};
end
end
