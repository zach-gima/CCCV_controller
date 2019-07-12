function plot_CCCV_data(out)
    fs = 25;
    
    current = out.yout{3}.Values;
    V = out.yout{1}.Values;
    
    figure('Position', [100 100 900 700])
    subplot(2,1,1)
    plot(current)
    ytickformat('%.2f')
    set(gca,'FontSize',fs)

    subplot(2,1,2)
    plot(V)
    ytickformat('%.2f')
    set(gca,'FontSize',fs)

end