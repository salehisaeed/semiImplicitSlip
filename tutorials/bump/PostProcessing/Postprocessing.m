% ============================ Initializing ============================= %
clc;
clear;
close all;
set(0, 'DefaultFigureRenderer', 'painters');
% set(0, 'DefaultFigureRenderer', 'OpenGL');
% ======================================================================= %



% ================================ Time  ================================ %
LineStyles = {'-','--','-.',':','-'};
Colors = [0 0 1; 1 0 0; 0 0 0; 0 0.5 0; 1 0 1];
lineWeight = [1 1 1 1.5 1.5];
picOutType = 'epsc';
Cases = {'slip','slipFixedNormal','slipSTL','slipSemiImplicitLaplacian'};
Fields = {'aspectRatio','nonOrthoAngle','cellAspectRatio','cellVolume','skewness'};
yLabs = {'Maximum cell aspect ratio','Maximum non-orthogonal angle'};
legNames = {'Slip','Fixed normal slip','Slip on surface','Semi-implicit slip'};
yLims = [1, 100; 20, 180];

postTime = [7 7 5.95 6.525];
negVolTime = [1.1, 2.775 2.25, 6.475];

for i = 1 : 2 % length(Fields)

    figure1 = figure;
    axes1 = axes('Parent',figure1,'FontSize',13,'YMinorTick','on',...
        'XMinorTick','on','TickLength',[0.02 0.005],'LineWidth',0.6);
    if i == 1
        set(axes1,'YScale','log');
    end

    box(axes1,'on');
    hold(axes1,'all');
    nSkip = 1;

    for j = 1 : length(Cases)

        Data = importdata(strcat('../Simulation/',Cases{j},...
            '/postProcessing/fieldMinMax_',Fields{i},'/',...
            num2str(postTime(j)),'/fieldMinMax.dat'));
        Data = Data.data;
        smoothData = [Data(:,1),smooth(Data(:,3),5)];

        if j < 4
            h(j) = plot(smoothData(:,1),smoothData(:,2),'LineStyle',...
                LineStyles{j},'LineWidth',1,'Color',Colors(j,:));
        elseif j < 5
            plot(smoothData(:,1),smoothData(:,2),'LineStyle','-',...
                'LineWidth',1,'Color',Colors(j,:));
            plot(smoothData(1:8:end,1),smoothData(1:8:end,2),'LineStyle',...
                'none','LineWidth',1,'Marker','.','Markersize',12,...
                'Color',Colors(j,:));
            h(j) = plot(inf,inf,'LineStyle','-','LineWidth',1,...
                'Marker','.','Markersize',12,'Color',Colors(j,:));
        else
            plot(smoothData(:,1),smoothData(:,2),'LineStyle','--',...
                'LineWidth',1,'Color',Colors(j,:));
            plot(smoothData(1:8:end,1),smoothData(1:8:end,2),'LineStyle',...
                'none','LineWidth',1,'Marker','^','Markersize',3.5,...
                'Color',Colors(j,:),'MarkerFaceColor',Colors(j,:));
            h(j) = plot(inf,inf,'LineStyle','--','LineWidth',1,...
                'Marker','^','Markersize',3.5,'Color',Colors(j,:));
        end

        plot([negVolTime(j) negVolTime(j)],[1,1e3],'LineStyle',':',...
            'LineWidth',1,'Color',Colors(j,:))
    end

    legend1 = legend(axes1,'show',h,legNames);
    set(legend1,'Location','Best','FontSize',13);
    xlim([0 6.4])
    ylim(yLims(i,:));
    xlabel('Time (s)','FontSize',15);
    ylabel(yLabs{i},'FontSize',15);

    saveas(figure1, Fields{i},'epsc')
end
% ======================================================================= %
