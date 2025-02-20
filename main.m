clear
close all

% Use full paths when possible to be robust
filename = mfilename('fullpath');
filepath = fileparts( filename );

addpath([filepath '/mesh2d']); initmsh();

% Medium-scale high-geometric-complexity aeropropulsive problem %%%%%%%%%%%%%%
surfaceFiles = {'dengwirda/mainElement.dat','onr-dep/nacelleVec30.dat', ...
    'dengwirda/foreFlap.dat','dengwirda/aftFlap.dat'};

for i = numel(surfaceFiles):-1:1
    fid = fopen([filepath '/airfoils/' surfaceFiles{i}],'r');
    surfaces{i} = cell2mat(textscan(fid,'%f%f','Delimiter',{'\t',','}));
    fclose(fid);
end

% Shift the nacelle to an appropriate position
surfaces{2} = surfaces{2} + [-0.1471 0.073046];


% Solve
opts.NumPanels = 200; % optionally pass options to the wake solver
[Cp,xc] = panel2d(surfaces,0,1,0.67,opts,'Plot','on','Colormap','vik');


% Plot Cp distributions
figure;
hold on;
for i = 1:numel(Cp)
    plot(xc{i},Cp{i})
end
set(gca,'YDir','reverse')