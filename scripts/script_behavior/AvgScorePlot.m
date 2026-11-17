%% prepare for sepsis scoring stair plot
%Wan-Ting Zhao, 2020-01-30 Universitaetsklinikum Jena
NormalDiet=zeros(24,1);%24 timepoint
FatDiet=zeros(24,1);
ProteinDiet=zeros(24,1);
NormalDiet= [0.00	2.50	3.00	3.40	3.00	3.20	3.40	3.80	4.00	4.20	3.60	3.60	3.60	1.80	1.00	0.80	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00];%avg scoring excluding MRT animals
FatDiet=[0.00	1.54	3.00	3.58	5.18	5.90	5.25	5.50	5.38	4.88	4.43	2.33	2.33	1.67	0.50	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00];
ProteinDiet=[0.00	2.00	4.09	3.38	3.75	4.75	4.00	2.40	2.40	2.60	2.60	1.80	1.80	1.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00];

NormalDietStd=zeros(24,1);%24 timepoint
FatDietStd=zeros(24,1);
ProteinDietStd=zeros(24,1);
NormalDietStd=[0.00	1.80	0.90	1.50	1.00	1.30	1.70	3.40	3.40	3.30	2.60	2.60	2.60	2.40	1.70	1.80	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00];
FatDietStd=[0.00	1.85	1.63	2.11	1.99	1.97	2.19	2.51	2.67	2.90	2.82	1.63	1.63	1.37	1.22	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00];
ProteinDietStd=[0.00	0.74	1.68	1.74	2.37	3.11	3.02	2.32	2.32	2.66	2.66	2.58	2.58	2.23	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00];


% %if sample is small use standard error instead of standard deviation.
% NormalDietStd=NormalDietStd./(length(NormalDietStd));
% FatDietStd=FatDietStd./(length(FatDietStd));
% ProteinDietStd=ProteinDietStd./(length(ProteinDietStd));

%set for xaxis value
xaxis(1:13)=0:0.25:3;
xaxis(14:24)=4:1:14;


%plot stair plot
figure1=figure;
axes1 = axes('Parent',figure1);
errorbar(xaxis,NormalDiet,NormalDietStd,'Color','g','LineWidth',3,'Marker','square','MarkerSize',12,'MarkerFaceColor','g');hold on;
errorbar(xaxis,FatDiet,FatDietStd,'Color','blue','LineWidth',3,'Marker','square','MarkerSize',12,'MarkerFaceColor','b');
errorbar(xaxis,ProteinDiet,ProteinDietStd,'Color','red','LineWidth',3,'Marker','square','MarkerSize',12,'MarkerFaceColor','r');


lgd = legend('Normal (n=9)','Fat (n=15)','Protein (n=12)');
lgd.FontSize = 15;
lgd.FontWeight = 'bold';

xlabel('Time span post infection (days)','FontSize',20)
ylabel('Disease (CSS)','FontSize',20)

axis([0 9 0 9])
 set(axes1,'FontSize',14,'FontWeight','bold','XTick',...
    [0 1 2 3 4 5 6 7 8 9 10 12 14],'YTick',...
    [0 1 2 3 4 5 6 7 8 9 10],'ZTick',[-1 0 1]);

grid on;

saveas(gcf, sprintf('AvgScorePlot%s.png',datestr(now,'yyyymmdd')));
saveas(gcf, sprintf('AvgScorePlot%s.fig',datestr(now,'yyyymmdd'))); 