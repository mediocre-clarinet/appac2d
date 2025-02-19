clear all
close all

addpath('./mesh2d'); initmsh();

surfaceFiles = {'mainVec0.dat','nacelleVec0.dat'};
for i = numel(surfaceFiles):-1:1
    fid = fopen(['airfoils/' surfaceFiles{i}],'r');
    surfaces{i} = cell2mat(textscan(fid,'%f%f','Delimiter',{'\t',','}));
    fclose(fid);
end

opts.FunctionTolerance = 1e-5;
[Cp,xc] = panel2d(surfaces,5,1,0.8,opts);

figure;
hold on;
for i = 1:numel(Cp)
    plot(xc{i},Cp{i})
end
set(gca,'YDir','reverse')