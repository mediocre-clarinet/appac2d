function flowVis(options,foils,varargin)
% FLOWVIS

% Meshing settings (not user defined)
opts.kind = 'delfront';
opts.rho2 = 1;


figure;
cmap = crameri(options.Colormap);
colormap(cmap);
hold on; axis image off;


gs = cmap*[0.299;0.587;0.114]; % grayscale conversion
LineColor = [0.8 0.8 0.8] - 0.6*round(gs(ceil(numel(gs)/2)));

if strcmpi(options.Mesh,'on')
    EdgeColor = LineColor;
else
    EdgeColor = 'none';
end


if nargin ~= 2
    wakes = varargin{1};
    k1 = varargin{2};
    k2 = varargin{3};

    N = wakes.m(1);
    % Distinguish jet and outer flow eval by making one wake left-running
    wakes.xo(1:N) = wakes.xo(1:N) + wakes.dx(1:N);
    wakes.yo(1:N) = wakes.yo(1:N) + wakes.dy(1:N);
    wakes.dx(1:N) = -wakes.dx(1:N);
    wakes.dy(1:N) = -wakes.dy(1:N);
    wakes.theta(1:N) = atan2(wakes.dy(1:N),wakes.dx(1:N));

    %--------------------------------------- Mesh jet volume
    node = [foils.co(k1+1:foils.m(1),:); ...
            wakes.co(1:N-1,:); ...
            flipud(wakes.co(N+1:2*N-1,:)); ...
            foils.co(foils.m(1)+(1:k2-1),:)];
%    node([false;node(1:end-1,1)>=4],:) = [];

    edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;

    [vert{1},~,tria{1},~] = refine2(node,edge,[],opts);

    [U,V] = influence(vert{1},foils,1);
    u{1} = U*foils.gamma + 1;
    v{1} = V*foils.gamma;
    [U,V] = influence(vert{1},wakes,-1); % -1 for jet interior
    u{1} = u{1} + U*wakes.gamma;
    v{1} = v{1} + V*wakes.gamma;

    %------------------------------------- Mesh outer volume
    node = [flipud(wakes.co(1:N-1,:)); ...
            foils.co(1:k1+1,:); ...
            foils.co(foils.m(1)+(k2-1:foils.m(2)),:); ...
            wakes.co(N+1:2*N-1,:); ...
            wakes.co(2*N-1,1) 2; ...
            -2 2;-2 -2; ...
            wakes.co(N-1,1) -2];

    edge = (1:size(node,1)).' + [0 1]; edge(end,2) = 1;

    if numel(foils.m) > 2
        j = cumsum(foils.m(3:end));
        edge = [edge;(1:j(end)).'+[0 1]+size(node,1)];
        edge(size(node,1)+j,2) = edge(size(node,1)+j,2) - foils.m(3:end).';
        node = [node;foils.co(foils.m(1)+foils.m(2)+1:end,:)];
    end

    [vert{2},~,tria{2},~] = refine2(node,edge,[],opts);

    [U,V] = influence(vert{2},foils,1);
    u{2} = U*foils.gamma + 1;
    v{2} = V*foils.gamma;
    [U,V] = influence(vert{2},wakes,1);
    u{2} = u{2} + U*wakes.gamma;
    v{2} = v{2} + V*wakes.gamma;

    %---------------------------------- Merge the two meshes
    TRI = [tria{1};tria{2}+size(vert{1},1)];
    VTX = [vert{1};vert{2}];
    data.u = [u{1};u{2}]; data.v = [v{1};v{2}];
else
    j = cumsum(foils.m);
    node = [foils.co;4 2;-2 2;-2 -2;4 -2];
    edge = (1:j(end)+4).' + [0 1]; edge(j,2) = edge(j,2) - foils.m.';
    edge(end,2) = j(end) + 1;
    [VTX,~,TRI,~] = refine2(node,edge,[],opts);
    [U,V] = influence(VTX,foils,1);
    data.u = U*foils.gamma + 1;
    data.v = V*foils.gamma;
end

data.q = sqrt(data.u.^2 + data.v.^2);
data.p = 1 - data.q.^2;
if nargin ~= 2
    CT = 0.5*((abs(wakes.gamma(N+1)) + 1)^2 - 1);
    data.p(1:numel(u{1})) = data.p(1:numel(u{1})) + 2*CT;
end

if strcmpi(options.CData,'p') || strcmpi(options.CData,'v')
    pivot = 0;
else
    pivot = 1;
end

% Create the contour plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patch('Faces',TRI,'Vertices',VTX, ...
    'FaceVertexCData',data.(lower(options.CData)), ...
    'FaceColor','interp','EdgeColor',EdgeColor);
cl = get(gca,'CLim');
set(gca,'CLim',max(abs(cl-pivot))*[-1 1]+pivot);

% Draw streamlines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FlowP = tristream(TRI,VTX(:,1),VTX(:,2),data.u,data.v,zeros(1,51)-1.99,-1.5:0.05:1);
for i = 1:numel(FlowP)
    plot(FlowP(i).x,FlowP(i).y,'-','Color',LineColor);
end