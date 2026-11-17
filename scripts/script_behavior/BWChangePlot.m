%% BW changes through sepsis progression plot
%Wan-Ting Zhao, 2020-01-30 Universitaetsklinikum Jena
NormalDiet=zeros(15,1);
FatDiet=zeros(15,1);
ProteinDiet=zeros(15,1);
xaxis=zeros(1,15);
NormalDiet= [0; -2;   -4;  -3;  -2;  -2;  0;   3;   5;   6;   7;  9;    9;  11;   12; ];
FatDiet=[0;  -1;  -4;  -5;  -3;  -2;  3;   8;   9;   10;  15;  16;  18;  18;  19;];
ProteinDiet=[0;  0;   -3;  -3;  -5;  -5;  -4;  4;   3;   2;   5;   6;   6;   7;   8;];

%set for xaxis value
xaxis(1:4)=[0  0.5 1.5 2.5];
xaxis(5:15)=4:1:14;

figure,plot(xaxis, NormalDiet,'b-o');hold on;
plot(xaxis, FatDiet,'r-o');
plot(xaxis, ProteinDiet,'k-o');
legend('Normal Diet','Fat Diet','Protein Diet');
legend('FontSize',15);

xlabel('Time span post infection (days)','FontSize',27)
ylabel('BW changes (%)','FontSize',27)


