clear all
close all

surfaceFiles = {'mainVec0.dat','nacelleVec0.dat'};
for i = numel(surfaceFiles):-1:1
    fid = fopen(surfaceFiles{i},'r');
    coords = fscanf(fid,'%f %f',[2 Inf]);
    surfaces{i} = coords.';
    fclose(fid);
end

[foils,wakes,Cp,xc] = appac2(surfaces,5,0.8,1,[]);

figure
plot(foils.xo,foils.yo,'k.')
hold on
plot(reshape(wakes.xo,[],2),reshape(wakes.yo,[],2),'b-')
set(gca,'DataAspectRatio',[1 1 1],'XLim',[-0.2 2])

figure
plot(xc,Cp)
set(gca,'YDir','reverse')

figure
g = reshape(wakes.gamma,[],2);
plot(reshape(wakes.xo,[],2),g(1:end-1,:))