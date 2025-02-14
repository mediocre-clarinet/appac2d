clear
close all

surfaceFiles = {'mainVec0.dat','nacelleVec0.dat'};
for i = numel(surfaceFiles):-1:1
    fid = fopen(surfaceFiles{i},'r');
    coords = fscanf(fid,'%f %f',[2 Inf]);
    surfaces{i} = coords.';
    fclose(fid);
end

[foils,surfaceVortexSheets,wakes,Cp,xc] = appac2(surfaces,10,0.8,1.0,[]);

figure
hold on
for i = 1:numel(foils)
    plot(foils(i).x,foils(i).y,'k-')
    plot(surfaceVortexSheets(i).x,surfaceVortexSheets(i).y,'r-','LineWidth',2)
    plot(wakes(i).x,wakes(i).y,'b-')
end
set(gca,'DataAspectRatio',[1 1 1],'XLim',[-0.2 2])

figure
plot(xc,Cp)
set(gca,'YDir','reverse')

%af = initSurface(coords.');
%
%[u,v] = influence(af.co,af);
%A = -u.*sin(af.theta) + v.*cos(af.theta);
%B = u.*cos(af.theta) + v.*sin(af.theta);
%A(af.m+1,1) = 1; A(end) = 1;
%RHS = [sin(af.theta);0];
%
%gamma = A \ RHS;
%
%Qtan = B*gamma + cos(af.theta);
%Cp = 1 - Qtan.^2;
%Cl = -Cp.'*af.dx;
%Cd = Cp.'*af.dy;
%
%figure
%plot(af.co(:,1),Cp)
%set(gca,'YDir','reverse')