function terminalSetPlot(font,X_set,Feasibleset_x0,INV_set,xk)
    figure;
    plot(X_set,'color','green','Alpha',0.5);
    hold on;
    plot(Feasibleset_x0,'Alpha',0.75);
    plot(INV_set,'color','yellow','Alpha',0.9);
    plot(xk(1,:),xk(2,:),'o')
    xlabel('$\bar{x}_1$ (voltage)',Interpreter="latex",Fontsize=font);
    ylabel('$\bar{x}_2$ (current)',Interpreter="latex",Fontsize=font);
    title('Feasible and terminal set',Interpreter="latex",Fontsize=font+4)
    legend({'State constraints','Feasible set', 'Terminal set'},Interpreter="latex",Fontsize=font)
end