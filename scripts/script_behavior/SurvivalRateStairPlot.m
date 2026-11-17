%% prepare for sepsis scoring stair plot
%Wan-Ting Zhao, 2020-01-30 Universitaetsklinikum Jena
NormalDiet=zeros(24,1);%24 timepoint
FatDiet=zeros(24,1);
ProteinDiet=zeros(24,1);
NormalDiet= [9	8	8	8	5	5	5	5	5	5	5	5	5	5	5	5	4	4	4	4	4	4	4	4];%survivor counts excluding MRT animals
FatDiet=[15	14	14	13	11	10	8	8	8	8	7	6	6	6	6	5	5	5	5	5	5	5	5	5];
ProteinDiet=[12	12	12	9	9	9	8	6	6	6	6	6	6	6	5	5	5	5	5	5	5	5	5	5];

%set for xaxis value
xaxis(1:13)=0:0.25:3;
xaxis(14:24)=4:1:14;

% prepare fraction of survival rate
FatDietFraction=FatDiet./FatDiet(1).*100;
NormalDietFraction=NormalDiet./NormalDiet(1).*100;
ProteinDietFraction=ProteinDiet./ProteinDiet(1).*100;

%plot stair plot
figure1=figure;
axes1 = axes('Parent',figure1);

stairs(xaxis,FatDietFraction,'Color','blue','LineWidth',3,'Marker','square');hold on;
stairs(xaxis,ProteinDietFraction,'Color','red','LineWidth',3,'Marker','square');
stairs(xaxis,NormalDietFraction,'Color','g','LineWidth',3,'Marker','square');

lgd = legend('Fat (n=15)','Protein (n=12)','Normal (n=9)');
lgd.FontSize = 15;
lgd.FontWeight = 'bold';

xlabel('Time span post infection (days)','FontSize',20)
ylabel('survival rate (%)','FontSize',20)

axis([0 14 0 100])
 set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
    [0 2 4 6 8 10 12 14],'YTick',...
    [0 20 40 60 80 100],'ZTick',[-1 0 1]);

grid on;

saveas(gcf, sprintf('SurvivalRateStairPlot%s.png',datestr(now,'yyyymmdd')));
% saveas(gcf, sprintf('SurvivalRateStairPlot%s.fig',datestr(now,'yyyymmdd'))); 