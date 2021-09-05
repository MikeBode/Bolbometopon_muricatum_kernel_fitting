function Plot_kernel(New)

if nargin == 0
    New = 1;
end

load OUTPUTS_random_size
% load OUTPUTS

xx = linspace(0,1000,200);
yy_bf = F(xx,Best_fit_k,Best_kernel);
CDF = cumsum(yy_bf./sum(yy_bf));
[~,Fifty] = min(abs(CDF-0.5));
[~,NinetyFive] = min(abs(CDF-0.95));
yy_bf = yy_bf./yy_bf(1);

if New == 1
    for i = 1:b
        yy(i,:) = F(xx,Bootstrap_k(i),Best_kernel);
        yy(i,:) = yy(i,:)./yy(i,1);
    end
    Q = quantile(yy,[0.01 0.95]);
    % Q(1,:) = F(xx,Confidence_bounds(1),Best_kernel); Q(1,:) = Q(1,:)./Q(1,1);
    % Q(2,:) = F(xx,Confidence_bounds(3),Best_kernel); Q(2,:) = Q(2,:)./Q(2,1);
    
    figure(1), clf, hold on, box off;
    BT = 0;
end
FS = 16; set(gcf,'color','w');
if New == 0
    BT = 9.15;
    LS = 155.3;
    xxx = xx;
    xx = xx./110;
    xx = xx+LS;
    yy_bf = yy_bf/2-BT;
    
    REM = find(xx>157);
    xx(REM) = [];
    yy_bf(REM) = [];
else
    BT = 0;
    LS = 0;
end

plot([LS LS + 2],[0 0]-BT,'k')
plot([LS LS],[0 0.6]-BT,'k')
plot(xx,yy_bf,'k','linewidth',1.5)
P = patch(xx([1 NinetyFive NinetyFive:-1:1]),[0-BT 0-BT yy_bf(NinetyFive:-1:1)],'k'); set(P,'facealpha',0.3,'edgecolor','none');
P = patch(xx([1 Fifty Fifty:-1:1]),[0-BT 0-BT yy_bf(Fifty:-1:1)],'k'); set(P,'facealpha',0.5,'edgecolor','none');
if New == 1
    plot(xx,Q(1,:),'k--')
    plot(xx,Q(2,:),'k--')
    xlabel('Distance (km)','fontsize',FS,'interpreter','latex')
    ylabel('Relative connectivity','fontsize',FS,'interpreter','latex')
else
	text(xx(Fifty),-BT-0.08,[num2str(round(xxx(Fifty))) 'km'],'fontsize',FS-7,'interpreter','latex','horizontalalignment','center')
	text(xx(NinetyFive),-BT-0.08,[num2str(round(xxx(NinetyFive))) 'km'],'fontsize',FS-7,'interpreter','latex','horizontalalignment','center')
    set(gca,'xticklabel',{},'yticklabel',{})
end

if New == 1
    set(gca,'fontsize',FS-2,'ytick',[0:0.2:1],'xtick',[0 round(xx(Fifty)) 50 round(xx(NinetyFive)) 100 150],'TickLabelInterpreter','latex')
    xlim([0 180])
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 25 12]*1.1)
    set(gcf, 'renderer', 'painters')
    print('-dtiff','-r200',['../../Figures/BestFitKernel.tiff'])
    print('-depsc',['../../Figures/BestFitKernel.eps'])
    print('-djpeg','-r600',['../../Figures/BestFitKernel.jpg'])
end


