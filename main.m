clear all
close all

surfaceFiles = {'mainVec0.dat','nacelleVec0.dat'};
for i = numel(surfaceFiles):-1:1
    fid = fopen(['airfoils/' surfaceFiles{i}],'r');
    coords = fscanf(fid,'%f %f',[2 Inf]);
    surfaces{i} = coords.';
    fclose(fid);
end

[foils,wakes,Cp,xc] = panel2d(surfaces,5,0.8,1);

figure
plot(xc,Cp)
set(gca,'YDir','reverse')