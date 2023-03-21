function trajectoryPlot(font,k_sim,xk,uk,xmin,xmax,umin,umax,shift)
    s = get(0, 'ScreenSize');
    figure('Position', [10 50 1200 500]);
    subplot(1,2,2);
    plot(0:k_sim,xk(:,1:length(xk)-1)',LineWidth=2);
    %hold on;
    %plot([0 k_sim],[xmax(1) xmax(1)],'color',"#0072BD",LineStyle="--",LineWidth=2);
    %plot([0 k_sim],[xmin(1) xmin(1)],'color',"#0072BD",LineStyle="--",LineWidth=2);
    %plot([0 k_sim],[xmax(2) xmax(2)],'color',"#D95319",LineStyle="--",LineWidth=2);
    %plot([0 k_sim],[xmin(2) xmin(2)],'color',"#D95319",LineStyle="--",LineWidth=2);
    xlim([0 k_sim]);
    xlabel('$k$',Interpreter="latex",Fontsize=font);
    if(strcmp('shifted',shift))
        ylabel('$\bar{x}_1$ [V], $\bar{x}_2$ [A]',Interpreter="latex",Fontsize=font);
        legend({'Voltage $\bar{x}_1$','Current $\bar{x}_2$'},Interpreter="latex",Fontsize=font);
    else
        ylim([xmin(1) xmax(1)])
        ylabel('$x_1$ [V], $x_2$ [A]',Interpreter="latex",Fontsize=font);
        legend({'Voltage $x_1$','Current $x_2$'},Interpreter="latex",Fontsize=font);
    end
    
    subplot(1,2,1);
    stairs(0:k_sim,uk,LineWidth=2);
    xlim([0 k_sim]);
    ylim([umin umax]);
    xlabel('$k$',Interpreter="latex",Fontsize=font);
    if(strcmp('shifted',shift))
        ylabel('duty cycle $\bar{u}$ [-]',Interpreter="latex",Fontsize=font);
    else
        ylabel('duty cycle $u$ [-]',Interpreter="latex",Fontsize=font);
    end
    if(strcmp('shifted',shift))
        sgtitle('Input and state trajectories or the shifted system',Interpreter="latex",Fontsize=font+4)
    else
        sgtitle('Input and state trajectories of the original non-shifted system',Interpreter="latex",Fontsize=font+4)
    end
end