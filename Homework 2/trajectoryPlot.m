function trajectoryPlot(font,k_sim,xk,uk,umin,umax,shift)
    figure;
    subplot(1,2,1);
    plot(0:k_sim,xk(:,1:length(xk)-1)',LineWidth=2);
    xlim([0 k_sim]);
    xlabel('$k$',Interpreter="latex",Fontsize=font);
    if(strcmp('shifted',shift))
        ylabel('$\bar{x}$',Interpreter="latex",Fontsize=font);
        legend({'$\bar{x}_1$ (voltage)','$\bar{x}_2$ (current)'},Interpreter="latex",Fontsize=font);
    else
        ylabel('$x$',Interpreter="latex",Fontsize=font);
        legend({'$x_1$ (voltage)','$x_2$ (current)'},Interpreter="latex",Fontsize=font);
    end
    
    subplot(1,2,2);
    stairs(0:k_sim,uk,LineWidth=2);
    xlim([0 k_sim]);
    ylim([umin-0.1 umax+0.1]);
    xlabel('$k$',Interpreter="latex",Fontsize=font);
    if(strcmp('shifted',shift))
        ylabel('$\bar{u}$  (duty cycle)',Interpreter="latex",Fontsize=font);
    else
        ylabel('$u$ (duty cycle)',Interpreter="latex",Fontsize=font);
    end
    sgtitle('Constrained MPC for the buck converter',Interpreter="latex",Fontsize=font+4)
end