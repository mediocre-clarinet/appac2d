clear all
close all

addpath('./mesh2d'); initmsh();

surfaceFiles = {'mainVec0.dat','nacelleVec0.dat','krueger.dat'};
for i = numel(surfaceFiles):-1:1
    fid = fopen(['airfoils/' surfaceFiles{i}],'r');
    surfaces{i} = cell2mat(textscan(fid,'%f%f','Delimiter',{'\t',','}));
    fclose(fid);
end

[foils,wakes,Cp,xc] = panel2d(surfaces,5,1,0.8);

figure
plot(xc,Cp)
set(gca,'YDir','reverse')