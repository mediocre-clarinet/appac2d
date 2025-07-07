function [wakes,gamma,iter,E] = solveWake(foils,Ainv,RHS,CT,opts)
% SOLVEWAKE  Solve for the global circulation solution and wake shape using
% Shollenberger's algorithm:
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

% Attach flat wake to the trailing edge of each propulsive element %%%%%%%%%%%
N = opts.NumPanels + 1; % add far-field panel to the panel count
wakes.m = [N N];
for i = 2:-1:1
    k = (i-1)*N+(1:N);
    if strcmpi(opts.NodeSpacing,'cosine')
        wakes.xo(k,:) = foils.xo(1+(i-1)*foils.m(1)) + ...
            opts.WakeLengthChords*(1-cos(linspace(0,pi/2,N))).';
    elseif strcmpi(opts.NodeSpacing,'uniform')
        wakes.xo(k,:) = foils.xo(1+(i-1)*foils.m(1)) + ...
            linspace(0,opts.WakeLengthChords,N).';
    end
    wakes.yo(k,:) = foils.yo(1+(i-1)*foils.m(1)) + zeros(N,1);
    wakes.dx(k,:) = [diff(wakes.xo(k)); 1e3];
end
wakes.dy = zeros(2*N,1);
wakes.theta = zeros(2*N,1);
wakes.co = [wakes.xo+wakes.dx/2 wakes.yo];

% Create initial guess for wake circulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gammaInf = sqrt(2*CT + 1) - 1;
wakes.gamma = repelem([gammaInf;-gammaInf],N+1);
gnew = wakes.gamma;

if strcmpi(opts.Display,'iter') || strcmpi(opts.Display,'final')
    figure;
    hold on; axis image; k = 0;
    for i = 1:numel(foils.m)
        plot(foils.xo(k+[1:foils.m(i) 1]),foils.yo(k+[1:foils.m(i) 1]),'k-');
        k = k + foils.m(i);
    end
    h(1) = plot(wakes.xo(1:N),wakes.yo(1:N),'b-');
    h(2) = plot(wakes.xo(N+1:2*N),wakes.yo(N+1:2*N),'r-');
end

iter = 0;
E = 1;
while (E > opts.FunctionTolerance) && (iter < opts.MaxIterations)
    iter = iter + 1;

    if strcmpi(opts.Display,'iter')
        set(h(1),'YData',wakes.yo(1:N));
        set(h(2),'YData',wakes.yo(N+1:2*N));
        drawnow;
    end

    % Solve airfoil circulation distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [U,V] = influence(foils.co,wakes,1);
    A = -U.*sin(foils.theta) + V.*cos(foils.theta);
    gamma = Ainv*(RHS - [A*wakes.gamma; ...
                         -wakes.gamma(1); ...
                         -wakes.gamma(N+2); ...
                         zeros(numel(foils.m)-2,1)]);

    % Calculate induced velocities on wake boundaries %%%%%%%%%%%%%%%%%%%%%%%%
    [U,V] = influence(wakes.co,foils,1);
    u = U*gamma + 1;
    v = V*gamma;
    [U,V] = influence(wakes.co,wakes,0);
    u = u + U*wakes.gamma;
    v = v + V*wakes.gamma;

    % Update wake shape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wakes.dy = v./u.*wakes.dx;
    wakes.dy([N 2*N]) = 0; % far-field panels remain flat
    wakes.yo(  2:N  ) = wakes.yo(1)   + cumsum(wakes.dy(  1:N-1  ));
    wakes.yo(N+2:2*N) = wakes.yo(N+1) + cumsum(wakes.dy(N+1:2*N-1));
    % Handle wakes crossing
    if wakes.yo(N) >= wakes.yo(2*N)
        yd = wakes.yo(2*N) - wakes.yo(N) - 0.025;
        % Rebuild YO of the lower wake directly rather than through CUMSUM
        wakes.yo(2:N) = wakes.yo(2:N) + yd* ...
            (wakes.xo(2:N)-wakes.xo(1))/opts.WakeLengthChords;
        wakes.dy(1:N-1) = diff(wakes.yo(1:N));
    end
    wakes.theta = atan2(wakes.dy,wakes.dx);
    wakes.co(:,2) = wakes.yo + wakes.dy/2;

    % Update wake circulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = CT./sqrt(u.*u + v.*v);
    g(N+1:2*N) = -g(N+1:2*N);
    % It is unclear if interpolating over x (faster) is less accurate than
    % interpolating over panel length (slower)
    gnew(  2:N    ) = g(  1:N-1  ) + diff(g(  1:N  )).* ...
        wakes.dx(  1:N-1  )./(wakes.dx(  1:N-1  )+wakes.dx(  2:N  ));
    gnew(N+3:2*N+1) = g(N+1:2*N-1) + diff(g(N+1:2*N)).* ...
        wakes.dx(N+1:2*N-1)./(wakes.dx(N+1:2*N-1)+wakes.dx(N+2:2*N));
    gnew(1)   = 2*g(1)   - gnew(2);
    gnew(N+2) = 2*g(N+1) - gnew(N+3);
    % Decay circulation to the far-field value
    k = find(wakes.xo(  1:N  ) >= wakes.xo(N)   - opts.GammaDecayChords, 1);
    gnew(    k:N    ) =  gammaInf + ( gammaInf-gnew(k)    )* ...
        (wakes.xo(  k:N  )-wakes.xo(N)  )/(wakes.xo(N)  -wakes.xo(k)  +eps);
    k = find(wakes.xo(N+1:2*N) >= wakes.xo(2*N) - opts.GammaDecayChords, 1);
    gnew(N+1+k:2*N+1) = -gammaInf + (-gammaInf-gnew(k+N+1))* ...
        (wakes.xo(N+k:2*N)-wakes.xo(2*N))/(wakes.xo(2*N)-wakes.xo(k+N)+eps);

    % Calculate residual between current and previous iteration %%%%%%%%%%%%%%
    E = sum(abs(gnew - wakes.gamma))/(2*N*gammaInf);

    % Adjust wake circulation by a relaxation factor %%%%%%%%%%%%%%%%%%%%%%%%%
    wakes.gamma = opts.RelaxationFactor*gnew + ...
        (1 - opts.RelaxationFactor)*wakes.gamma;
end

if strcmpi(opts.Display,'iter') || strcmpi(opts.Display,'final')
    set(h(1),'YData',wakes.yo(1:N));
    set(h(2),'YData',wakes.yo(N+1:2*N));
    drawnow;
end
