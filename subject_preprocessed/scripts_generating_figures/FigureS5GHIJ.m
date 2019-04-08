C_1A1=[180 176 196]/255;
C_1A11=[145 173 195]/255;
C_2=[241 220 212]/255;

%% 
figure(22)
subplot(2,3,1);hold on;
As=[];Bs=[];
for i = 1;
    SL1 = Ratios{i}.*ODS{i};
    As = [As SL1];
    SL2 = (1-Ratios{i}).*ODS{i};
    Bs = [Bs SL2];
end
SL1 = mean(As,2);
SL2 = mean(Bs,2);
h = area(tmQ,[SL1,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylim([0 0.35]);
h(1).FaceColor = C_1A11;
h(2).FaceColor = C_1A1;
xlabel('hours');
ylabel('Abunance (OD*Fraction)');
%%
subplot(2,3,4);hold on;
errorbar(tmQ, Ratio_1_mean_ts,Ratio_1_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_1_mean_ts,Ratio_1_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);

ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylabel('Fraction');xlabel('hours');

%% 
subplot(2,3,2);hold on;
As=[];Bs=[];
for i = 2;
    SL1 = Ratios{i}.*ODS{i};
    As = [As SL1];
    SL2 = (1-Ratios{i}).*ODS{i};
    Bs = [Bs SL2];
end
SL1 = mean(As,2);
SL2 = mean(Bs,2);
h = area(tmQ,[SL1,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylim([0 0.35]);
h(1).FaceColor = C_1A11;
h(2).FaceColor = C_1A1;
xlabel('hours');
ylabel('Abunance (OD*Fraction)');
%%
subplot(2,3,5);hold on;
errorbar(tmQ, Ratio_2_mean_ts,Ratio_1_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_2_mean_ts,Ratio_1_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);

ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylabel('Fraction');xlabel('hours');

%% 
subplot(2,3,3);hold on;
As=[];Bs=[];
for i = 3;
    SL1 = Ratios{i}.*ODS{i};
    As = [As SL1];
    SL2 = (1-Ratios{i}).*ODS{i};
    Bs = [Bs SL2];
end
SL1 = mean(As,2);
SL2 = mean(Bs,2);
h = area(tmQ,[SL1,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylim([0 0.35]);
h(1).FaceColor = C_1A11;
h(2).FaceColor = C_1A1;
xlabel('hours');
ylabel('Abunance (OD*Fraction)');
%%
subplot(2,3,6);hold on;
errorbar(tmQ, Ratio_3_mean_ts,Ratio_1_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_3_mean_ts,Ratio_1_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);

ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylabel('Fraction');xlabel('hours');
%%
for i = 1:6;
    subplot(2,3,i);
    set(gca,'FontSize',15);
    set(gca,'LineWidth',1);
    xlabel([])
    ylabel([])

end


%%
C_1A1=[180 176 196]/255;
C_1A11=[145 173 195]/255;
C_2=[241 220 212]/255;

%% 
figure(233445)
subplot(2,3,1);hold on;
As=[];Bs=[];
for i = 4:6:15;
    SL1 = Ratios{i}.*ODS{i};
    As = [As SL1];
    SL2 = (1-Ratios{i}).*ODS{i};
    Bs = [Bs SL2];
end
SL1 = mean(As,2);
SL2 = mean(Bs,2);
h = area(tmQ,[SL1,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylim([0 0.35]);
h(1).FaceColor = C_1A1;
h(2).FaceColor = C_2;
xlabel('hours');
ylabel('Abunance (OD*Fraction)');
%%
subplot(2,3,4);hold on;
errorbar(tmQ, Ratio_4_mean_ts,Ratio_4_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ, Ratio_10_mean_ts,Ratio_10_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_4_mean_ts,Ratio_4_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_10_mean_ts,Ratio_10_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylabel('Fraction');xlabel('hours');

%%
%% 
subplot(2,3,2);hold on;
ylim([0 0.35]);

As=[];Bs=[];
for i = 5:6:15;
    SL1 = Ratios{i}.*ODS{i};
    As = [As SL1];
    SL2 = (1-Ratios{i}).*ODS{i};
    Bs = [Bs SL2];
end
SL1 = mean(As,2);
SL2 = mean(Bs,2);
Tot = SL1+SL2;
h = area(tmQ,[SL1,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
h(1).FaceColor = C_1A1;
h(2).FaceColor = C_2;
xlabel('hours');
ylabel('Abunance (OD*Fraction)');
%%
subplot(2,3,5);hold on;
errorbar(tmQ, Ratio_5_mean_ts,Ratio_5_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ, Ratio_11_mean_ts,Ratio_11_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_5_mean_ts,Ratio_5_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_11_mean_ts,Ratio_11_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
xlabel('hours');
ylabel('Fraction');
%% 
subplot(2,3,3);hold on;
As=[];Bs=[];
ylim([0 0.35]);

for i = 6:6:15;
    SL1 = Ratios{i}.*ODS{i};
    As = [As SL1];
    SL2 = (1-Ratios{i}).*ODS{i};
    Bs = [Bs SL2];
end
SL1 = mean(As,2);
SL2 = mean(Bs,2);
Tot = SL1+SL2;
h = area(tmQ,[SL1,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
h(1).FaceColor = C_1A1;
h(2).FaceColor = C_2;
xlabel('hours');
ylabel('Abunance (OD*Fraction)');

%%
subplot(2,3,6);hold on;
errorbar(tmQ, Ratio_6_mean_ts,Ratio_6_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ, Ratio_12_mean_ts,Ratio_12_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_6_mean_ts,Ratio_6_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_12_mean_ts,Ratio_12_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
xlabel('hours');
ylabel('Fraction');
%%
for i = 1:6;
    subplot(2,3,i);
    set(gca,'FontSize',15);
    set(gca,'LineWidth',1);
    xlabel([])
    ylabel([])

end

% C_1A11=[180 176 196]/255;
C_1A11=[145 173 195]/255;
C_2=[241 220 212]/255;

%% 
figure(23)
subplot(2,3,1);hold on;
As=[];Bs=[];
for i = 7:6:15;
    SL1 = Ratios{i}.*ODS{i};
    As = [As SL1];
    SL2 = (1-Ratios{i}).*ODS{i};
    Bs = [Bs SL2];
end
SL1 = mean(As,2);
SL2 = mean(Bs,2);
h = area(tmQ,[SL1,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylim([0 0.35]);
h(1).FaceColor = C_1A11;
h(2).FaceColor = C_2;
xlabel('hours');
ylabel('Abunance (OD*Fraction)');
%%
subplot(2,3,4);hold on;
errorbar(tmQ, Ratio_7_mean_ts,Ratio_7_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ, Ratio_13_mean_ts,Ratio_13_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_7_mean_ts,Ratio_7_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_13_mean_ts,Ratio_13_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylabel('Fraction');xlabel('hours');

%%
%% 
subplot(2,3,2);hold on;
ylim([0 0.35]);

As=[];Bs=[];
for i = 8:6:15;
    SL1 = Ratios{i}.*ODS{i};
    As = [As SL1];
    SL2 = (1-Ratios{i}).*ODS{i};
    Bs = [Bs SL2];
end
SL1 = mean(As,2);
SL2 = mean(Bs,2);
Tot = SL1+SL2;
h = area(tmQ,[SL1,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
h(1).FaceColor = C_1A11;
h(2).FaceColor = C_2;
xlabel('hours');
ylabel('Abunance (OD*Fraction)');
%%
subplot(2,3,5);hold on;
errorbar(tmQ, Ratio_8_mean_ts,Ratio_8_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ, Ratio_14_mean_ts,Ratio_14_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_8_mean_ts,Ratio_8_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_14_mean_ts,Ratio_14_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
xlabel('hours');
ylabel('Fraction');
%% 
subplot(2,3,3);hold on;
As=[];Bs=[];
ylim([0 0.35]);

for i = 9:6:15;
    SL1 = Ratios{i}.*ODS{i};
    As = [As SL1];
    SL2 = (1-Ratios{i}).*ODS{i};
    Bs = [Bs SL2];
end
SL1 = mean(As,2);
SL2 = mean(Bs,2);
Tot = SL1+SL2;
h = area(tmQ,[SL1,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
h(1).FaceColor = C_1A11;
h(2).FaceColor = C_2;
xlabel('hours');
ylabel('Abunance (OD*Fraction)');

%%
subplot(2,3,6);hold on;
errorbar(tmQ, Ratio_9_mean_ts,Ratio_9_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ, Ratio_15_mean_ts,Ratio_15_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_9_mean_ts,Ratio_9_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
errorbar(tmQ, 1-Ratio_15_mean_ts,Ratio_15_std_ts,'Color',C_2,'LineWidth',2,'CapSize',10);
ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
xlabel('hours');
ylabel('Fraction');
%%
for i = 1:6;
    subplot(2,3,i);
    set(gca,'FontSize',15);
    set(gca,'LineWidth',1);
    xlabel([])
    ylabel([])

end

%%
C_1A11=[180 176 196]/255;
C_1A1=[145 173 195]/255;
C_2=[241 220 212]/255;

%% 
figure(234)
subplot(2,3,1);hold on;
As=[];Bs=[];
SL1a1x = (1-Ratio_16_1_mean_ts).*ODS{16+1};
SL1a1y = (1-Ratio_19_1_mean_ts).*ODS{19+1};
SL1a11x = Ratio_16_3_mean_ts.*ODS{16+1};
SL1a11y = Ratio_19_3_mean_ts.*ODS{19+1};
SL2x = (Ratio_16_2_mean_ts).*ODS{16+1};
SL2y = (Ratio_19_2_mean_ts).*ODS{19+1};


As = [SL1a1x SL1a1y];
Bs = [SL1a11x SL1a11y];
Cs = [SL2x SL2y];

SL1a1 = mean(As,2);
SL1a11 = mean(Bs,2);
SL2 = mean(Cs,2);
h = area(tmQ+randn(6,1)*0.09,[SL1a1,SL1a11,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylim([0 0.35]);
h(1).FaceColor = C_1A1;
h(2).FaceColor = C_1A11;
h(3).FaceColor = C_2;

xlabel('hours');
ylabel('Abunance (OD*Fraction within SL1)');
%%
subplot(2,3,4);hold on;
errorbar(tmQ+randn(6,1)*0.09, 1-Ratio_16_4_mean_ts,Ratio_16_4_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ+randn(6,1)*0.09, 1-Ratio_19_4_mean_ts,Ratio_19_4_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ+randn(6,1)*0.09, Ratio_16_4_mean_ts,Ratio_16_4_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ+randn(6,1)*0.09, Ratio_19_4_mean_ts,Ratio_19_4_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);

ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylabel('Fraction within SL1');xlabel('hours');


%% 
subplot(2,3,2);hold on;
As=[];Bs=[];
SL1a1x = (1-Ratio_17_1_mean_ts).*ODS{17+1};
SL1a1y = (1-Ratio_20_1_mean_ts).*ODS{20+1};
SL1a11x = Ratio_17_3_mean_ts.*ODS{17+1};
SL1a11y = Ratio_20_3_mean_ts.*ODS{20+1};
SL2x = (Ratio_17_2_mean_ts).*ODS{17+1};
SL2y = (Ratio_20_2_mean_ts).*ODS{20+1};


As = [SL1a1x SL1a1y];
Bs = [SL1a11x SL1a11y];
Cs = [SL2x SL2y];

SL1a1 = mean(As,2);
SL1a11 = mean(Bs,2);
SL2 = mean(Cs,2);
h = area(tmQ+randn(6,1)*0.09,[SL1a1,SL1a11,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylim([0 0.35]);
h(1).FaceColor = C_1A1;
h(2).FaceColor = C_1A11;
h(3).FaceColor = C_2;

xlabel('hours');
ylabel('Abunance (OD*Fraction within SL1)');
%%
subplot(2,3,5);hold on;
errorbar(tmQ+randn(6,1)*0.09, 1-Ratio_17_4_mean_ts,Ratio_17_4_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ+randn(6,1)*0.09, 1-Ratio_20_4_mean_ts,Ratio_20_4_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ+randn(6,1)*0.09, Ratio_17_4_mean_ts,Ratio_17_4_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ+randn(6,1)*0.09, Ratio_20_4_mean_ts,Ratio_20_4_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylabel('Fraction within SL1');xlabel('hours');

%% 
subplot(2,3,3);hold on;
As=[];Bs=[];
SL1a1x = (1-Ratio_18_1_mean_ts).*ODS{18+1};
SL1a1y = (1-Ratio_21_1_mean_ts).*ODS{21+1};
SL1a11x = Ratio_18_3_mean_ts.*ODS{18+1};
SL1a11y = Ratio_21_3_mean_ts.*ODS{21+1};
SL2x = (Ratio_18_2_mean_ts).*ODS{18+1};
SL2y = (Ratio_21_2_mean_ts).*ODS{21+1};


As = [SL1a1x SL1a1y];
Bs = [SL1a11x SL1a11y];
Cs = [SL2x SL2y];

SL1a1 = mean(As,2);
SL1a11 = mean(Bs,2);
SL2 = mean(Cs,2);
h = area(tmQ+randn(6,1)*0.09,[SL1a1,SL1a11,SL2]);
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylim([0 0.35]);
h(1).FaceColor = C_1A1;
h(2).FaceColor = C_1A11;
h(3).FaceColor = C_2;

xlabel('hours');
ylabel('Abunance (OD*Fraction within SL1)');
%%
subplot(2,3,6);hold on;
errorbar(tmQ+randn(6,1)*0.09, 1-Ratio_18_4_mean_ts,Ratio_18_4_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ+randn(6,1)*0.09, 1-Ratio_21_4_mean_ts,Ratio_21_4_std_ts,'Color',C_1A1,'LineWidth',2,'CapSize',10);
errorbar(tmQ+randn(6,1)*0.09, Ratio_18_4_mean_ts,Ratio_18_4_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
errorbar(tmQ+randn(6,1)*0.09, Ratio_21_4_mean_ts,Ratio_21_4_std_ts,'Color',C_1A11,'LineWidth',2,'CapSize',10);
ylim([0 1])
box on; set(gca,'TickDir','out');
xlim([0 22]);
ylabel('Fraction within SL1');xlabel('hours');

%%
%%
for i = 1:6;
    subplot(2,3,i);
    set(gca,'FontSize',15);
    set(gca,'LineWidth',1);
    xlabel([])
    ylabel([])

end
