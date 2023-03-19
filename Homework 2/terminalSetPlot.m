function terminalSetPlot(font,Feasibleset_x0,INV_set)
    figure;
    plot(Feasibleset_x0);
    hold on;
    plot(INV_set,'color','yellow');
    xlabel('$\bar{x}_1$ (voltage)',Interpreter="latex",Fontsize=font);
    ylabel('$\bar{x}_2$ (current)',Interpreter="latex",Fontsize=font);
    title('Feasible and terminal set',Interpreter="latex",Fontsize=font+4)
    legend({'Feasible set', 'Terminal set'},Interpreter="latex",Fontsize=font)
end