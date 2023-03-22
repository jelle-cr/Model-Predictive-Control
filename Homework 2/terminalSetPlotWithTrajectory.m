function terminalSetPlotWithTrajectory(font,X_set,Feasibleset_x0,INV_set,x1,x2,x3,x4)
    s = get(0, 'ScreenSize');
    figure('Position', [10 50 700 500]);
    plot(X_set,'color','green','Alpha',0.5);
    hold on;
    plot(Feasibleset_x0,'color','red');
    plot(INV_set,'color','yellow');
    plot(x3(1,:),x3(2,:),'Color','blue',LineWidth=2);
    plot(x2(1,:),x2(2,:),'Color','#016615',LineWidth=2);
    plot(x4(1,:),x4(2,:),'Color','#02bdde',LineWidth=2);
    plot(x1(1,:),x1(2,:),'Color','#fcb103',LineWidth=2);
    plot(x3(1,1),x3(2,1),'square','Color','blue',MarkerSize=8,LineWidth=2);
    plot(x3(1,length(x3)),x3(2,length(x3)),'^','Color','blue',MarkerSize=8,LineWidth=2);
    plot(x2(1,1),x2(2,1),'square','Color','#016615',MarkerSize=8,LineWidth=2);
    plot(x2(1,length(x2)),x2(2,length(x2)),'^','Color','#016615',MarkerSize=8,LineWidth=2);
    plot(x4(1,1),x4(2,1),'square','Color','#02bdde',MarkerSize=8,LineWidth=2);
    plot(x4(1,length(x4)),x4(2,length(x4)),'^','Color','#02bdde',MarkerSize=8,LineWidth=2);
    plot(x1(1,1),x1(2,1),'square','Color','#fcb103',MarkerSize=8,LineWidth=2);
    plot(x1(1,length(x1)),x1(2,length(x1)),'^','Color','#fcb103',MarkerSize=8,LineWidth=2);
    ax = gca;
    ax.FontSize=font-6;
    xlabel('$\bar{x}_1$ [V]',Interpreter="latex",Fontsize=font);
    ylabel('$\bar{x}_2$ [A]',Interpreter="latex",Fontsize=font);
    title('Feasible and terminal set with state trajectories',Interpreter="latex",Fontsize=font+2)
    %legend({'State constraints','Feasible set', 'Terminal set'},Interpreter="latex",Fontsize=font-2)
end