function terminalSetPlot(font,X_set,Feasibleset_x0,INV_set,xk)
    s = get(0, 'ScreenSize');
    figure('Position', [10 50 700 500]);
    plot(X_set,'color','green','Alpha',0.5);
    hold on;
    plot(Feasibleset_x0);
    plot(INV_set,'color','yellow');
    %plot(xk(1,:),xk(2,:),'o');
    ax = gca;
    ax.FontSize=font-6;
    xlabel('$\bar{x}_1$ [V]',Interpreter="latex",Fontsize=font);
    ylabel('$\bar{x}_2$ [A]',Interpreter="latex",Fontsize=font);
    title('Feasible and terminal set',Interpreter="latex",Fontsize=font+2)
    legend({'State constraints','Feasible set', 'Terminal set'},Interpreter="latex",Fontsize=font-2)
end