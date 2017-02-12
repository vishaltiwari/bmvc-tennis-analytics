function [] = TrajectoryPlot3D(XYZts)
    Dx=23.7744;Dy=9.6;Dn=0.914;Ds=5.4864;

    % lines
    hold on
    plot3(Dx*[0 1 1 0 0],Dy*[0 0 1 1 0],Dn*[0 0 0 0 0],'r'); % court
    plot3(0.5*Dx*[1 1 1 1 1],Dy*[0 0 1 1 0],Dn*[0 1 1 0 0],'r'); % net
    plot3(Ds*[1 1],Dy*[0 1],Dn*[0 0],'r'); % service line 1
    plot3((Dx-Ds)*[1 1],Dy*[0 1],Dn*[0 0],'r'); % service line 2
    plot3([Ds Dx-Ds],0.5*Dy*[1 1],Dn*[0 0],'r'); % half line
    hold off

    xlabel('x (m), ');
    ylabel('y (m), ');
    zlabel('h (m), ');
    title('tennis ball trajectory from 0 to ');
    
    hold on
    
    plot3(XYZts(:,1) , XYZts(:,2), XYZts(:,3), 'ko-', 'linewidth', 0.5)

    rotate3d on;
end