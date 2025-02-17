function [wakes,gamma,iter,E] = solveWake(foils,Ainv,RHS,CT,options)
% Calculate the global circulation solution using Shollenberger's algorithm:
%   1. Guess initial wake shape and bound circulation
%   2. Calculate induced velocities on the airfoil surfaces due to the wake
%   3. Solve circulation strengths on the airfoil surfaces such that there is
%      no flow through the solid boundaries
%   4. Calculate induced velocities on the wake
%   5. Move the wake panels to satisfy flow tangency on the wake boundaries
%   6. Calculate bound circulation on the wake from the newly calculated
%      average velocity across the wake boundary
%   7. If the new wake circulation is close to the previous solution, exit the
%      program, else return to step 2

% Default options
opts.MaxIterations = 50;
opts.FunctionTolerance = 1e-6;
opts.RelaxationFactor = 0.5;
opts.NumPanels = 200;
opts.WakeLengthChords = 7;
opts.GammaDecayChords = 2;
opts.NodeSpacing = 'cosine';
opts.Display = 'iter';

% Attach flat wake to the trailing edge of each propulsive element %%%%%%%%%%%
N = opts.NumPanels + 1; % add far-field panel to the panel count
wakes.m = [N N];
for i = 2:-1:1
    k = (i-1)*N+(1:N); % working indices
    if strcmpi(opts.NodeSpacing,'cosine')
        wakes.xo(k,:) = foils.xo((i-1)*foils.m(1)+1) + ...
            opts.WakeLengthChords*(1-cos(linspace(0,pi/2,N))).';
    elseif strcmpi(opts.NodeSpacing,'uniform')
        wakes.xo(k,:) = foils.xo((i-1)*foils.m(1)+1) + ...
            linspace(0,opts.WakeLengthChords,N).';
    end
    wakes.yo(k,:) = foils.yo((i-1)*foils.m(1)+1) + zeros(N,1);
    wakes.dx(k,:) = [diff(wakes.xo(k)); 1e3];
    wakes.dy(k,:) = zeros(N,1);
end
wakes.theta = atan2(wakes.dy,wakes.dx);
wakes.co = [wakes.xo+wakes.dx/2 wakes.yo+wakes.dy/2];

% Create initial guess for wake circulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammaInf = sqrt(2*CT + 1) - 1;
wakes.gamma = repelem([gammaInf;-gammaInf],N+1);

iter = 0;
E = 1;
while (E > opts.FunctionTolerance) && (iter < opts.MaxIterations)
    iter = iter + 1;

    % Solve airfoil circulation distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [U,V] = influence(foils.co,wakes,pi);
    A = -U.*sin(foils.theta) + V.*cos(foils.theta);
    gamma = Ainv*(RHS - [A*wakes.gamma; ...
                         -wakes.gamma(1); ...
                         -wakes.gamma(N+2); ...
                         zeros(numel(foils.m)-2,1)]);

    % Calculate induced velocities on wake boundaries %%%%%%%%%%%%%%%%%%%%%%%%
    [U,V] = influence(wakes.co,foils,pi);
    u = U*gamma + 1;
    v = V*gamma;
    [U,V] = influence(wakes.co,wakes,0); % 0 yields the average of pi and -pi
    u = u + U*wakes.gamma;
    v = v + V*wakes.gamma;

    % Update wake shape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wakes.dy = v./u.*wakes.dx;
    wakes.dy([N 2*N]) = 0; % far-field panels remain flat
    wakes.yo(  2:N  ) = wakes.yo(1)   + cumsum(wakes.dy(  1:N-1  ));
    wakes.yo(N+2:2*N) = wakes.yo(N+1) + cumsum(wakes.dy(N+1:2*N-1));
    wakes.theta = atan2(wakes.dy,wakes.dx);
    wakes.co(:,2) = wakes.yo + wakes.dy/2;

    % Update wake circulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gprev = wakes.gamma;
    g = CT./sqrt(u.^2 + v.^2);
    g(N+1:2*N) = -g(N+1:2*N);
    % It is unclear if interpolating over x (faster) is less accurate than
    % interpolating over panel length (slower)
    wakes.gamma(  2:N    ) = g(  1:N-1  ) + diff(g(  1:N  )).* ...
        wakes.dx(  1:N-1  )./(wakes.dx(  1:N-1  )+wakes.dx(  2:N  ));
    wakes.gamma(N+3:2*N+1) = g(N+1:2*N-1) + diff(g(N+1:2*N)).* ...
        wakes.dx(N+1:2*N-1)./(wakes.dx(N+1:2*N-1)+wakes.dx(N+2:2*N));
    wakes.gamma(1)   = 2*g(1)   - wakes.gamma(2);
    wakes.gamma(N+2) = 2*g(N+1) - wakes.gamma(N+3);
    % Decay circulation to the far-field value
    k = find(wakes.xo(1:N) >= wakes.xo(N) - opts.GammaDecayChords, 1);
    wakes.gamma(    k:N    ) =  gammaInf + ( gammaInf-wakes.gamma(k)    )* ...
        (wakes.xo(  k:N  )-wakes.xo(N)  )/(wakes.xo(N)  -wakes.xo(k)  +eps);
    k = find(wakes.xo(N+1:2*N) >= wakes.xo(2*N) - opts.GammaDecayChords, 1);
    wakes.gamma(k+N+1:2*N+1) = -gammaInf + (-gammaInf-wakes.gamma(k+N+1))* ...
        (wakes.xo(k+N:2*N)-wakes.xo(2*N))/(wakes.xo(2*N)-wakes.xo(k+N)+eps);

    % Calculate residual between current and previous iteration %%%%%%%%%%%%%%
    E = sum(abs(wakes.gamma - gprev))/(2*N*gammaInf);

    % Adjust wake circulation forwarded to next iteration by relaxation factor
    wakes.gamma = opts.RelaxationFactor*wakes.gamma + ...
        (1 - opts.RelaxationFactor)*gprev;
end
