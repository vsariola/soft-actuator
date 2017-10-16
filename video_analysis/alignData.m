function [alignedPressure,alignedCurvature] = alignData(pressure,curvature)
    plot(curvature);
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    title('Curvature. Click on the time when the experiment starts');
    leftc = ginput(1);
    title('Curvature. Click on the time when the experiment ends');
    rightc = ginput(1);
    plot(pressure);
    title('Pressure. Click on the time when the experiment starts');
    leftp = ginput(1);
    title('Pressure. Click on the time when the experiment ends');
    rightp = ginput(1);
    close(gcf);
    c1 = min(max(round(leftc(1)),1),length(curvature));
    c2 = min(max(round(rightc(1)),1),length(curvature));
    p1 = min(max(round(leftp(1)),1),length(pressure));
    p2 = min(max(round(rightp(1)),1),length(pressure));
    alignedCurvature = curvature(c1:c2);
    alignedPressure = pressure(round(linspace(p1,p2,c2-c1+1)));
end
    