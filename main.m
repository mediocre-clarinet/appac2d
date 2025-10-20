clear
close all

% Use full paths when possible to be robust
filename = mfilename('fullpath');
filepath = fileparts( filename );

addpath([filepath '/mesh2d']); initmsh();

% Medium-scale high-geometric-complexity aeropropulsive problem %%%%%%%%%%%%%%
surfaceFiles = {'dengwirda/mainElement.dat','dengwirda/nacelle.dat','dengwirda/aftFlap.dat'};

for i = numel(surfaceFiles):-1:1
    fid = fopen([filepath '/airfoils/' surfaceFiles{i}],'r');
    surfaces{i} = cell2mat(textscan(fid,'%f%f','Delimiter',{'\t',','}));
    fclose(fid);
end


% Solve
opts.NumPanels = 200; % optionally pass options to the wake solver
[Cp,xc] = panel2d(surfaces,5,1,0.67,opts,'Plot','on','Colormap','vik');


% Plot Cp distributions
figure;
hold on;
for i = 1:numel(Cp)
    plot(xc{i},Cp{i})
    xlabel('x/C')
    ylabel('Cp')

end
set(gca,'YDir','reverse')