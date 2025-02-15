clear all
close all

surfaceFiles = {'mainVec0.dat','nacelleVec0.dat'};
for i = numel(surfaceFiles):-1:1
    fid = fopen(surfaceFiles{i},'r');
    coords = fscanf(fid,'%f %f',[2 Inf]);
    surfaces{i} = coords.';
    fclose(fid);
end

[foils,surfaceVortexSheets,Cp,xc] = appac2(surfaces,10,0.8,1.5,[]);

figure
plot(foils.xo,foils.yo,'k.')
hold on
for i = 1:numel(surfaceVortexSheets)
    plot(surfaceVortexSheets(i).x,surfaceVortexSheets(i).y,'r-','LineWidth',2)
end
set(gca,'DataAspectRatio',[1 1 1],'XLim',[-0.2 2])

figure
plot(xc,Cp)
set(gca,'YDir','reverse')