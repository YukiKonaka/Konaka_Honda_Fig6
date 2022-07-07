clear all 
close all
clc

rng default
rng(3)

numberofparticles=100000;

% load animal data
C_noise=0.5;
load ('animal_data.mat')
N_step=length(a);


% run particle filter
myPF = particleFilter(@state_15_self_org,@like_15_self_org);

am0=1;
vw0 = 0.4;
Co = -5;
POo=0.5;
init=[0;0;0.1;0.1;Co;am0;vw0;POo];

F1 = init(1); 
F2 = init(2);
F3 = init(3); 
F4 = init(4);
F5 = init(5);
F6 = init(6);
F7 = init(7);
F8 = init(8);

F = [F1,F2,F3,F4,F5,F6,F7,F8];

S1 = (0.01)^2;
S2 = (0.01)^2;
S3 = (0.1)^2;
S4 = (0.1)^2;
S5 = (0.001);
S6 = (0.5);
S7 = (0.5);
S8 = (0.5);
S= diag([S1,S2,S3,S4,S5,S6,S7,S8]);

initialize(myPF,numberofparticles,F,S);
myPF.Particles(3:4,:)=gamrnd(10,0.001,[2,numberofparticles]);
myPF.Particles(5,:)=-15+30*rand(1,numberofparticles);
myPF.Particles(6,:)=0.04+0.02*rand(1,numberofparticles);
myPF.Particles(7,:)=0.2+0.5*rand(1,numberofparticles);
myPF.Particles(8,:)=0.3+0.4*rand(1,numberofparticles);


myPF.StateEstimationMethod = 'mean';
myPF.ResamplingMethod = 'systematic';

zEst=zeros(length(a),length(F));


for k=1:length(a)

    zEst(k,:) = correct(myPF,a(k),o(k));
 
    predict(myPF,a(k),o(k),C_noise);
    
    vv(:,:,k)=cov(transpose(myPF.Particles(1:5,:)));
        
end

v(:,:,:)=vv;



PF_ml=zEst(:,1);
PF_mr=zEst(:,2);
PF_pl=zEst(:,3);
PF_pr=zEst(:,4);
PF_C =zEst(:,5);
PF_am =zEst(:,6);
PF_vw =zEst(:,7);
PF_Po =zEst(:,8);

pf_ml=PF_ml;
pf_mr=PF_mr;
pf_pl=PF_pl;
pf_pr=PF_pr;
pf_C =PF_C;
pf_am =PF_am;
pf_vw =PF_vw;
pf_Po =PF_Po;

pf_v =v;
 
pf_rpl=1./(1+exp(-pf_ml));
pf_rpr=1./(1+exp(-pf_mr));
pf_precl=pf_pl./(pf_rpl.^2)./((1-pf_rpl).^2);
pf_precr=pf_pr./(pf_rpr.^2)./((1-pf_rpr).^2);

filter=[{pf_ml},{pf_mr},{pf_pl},{pf_pr},{pf_C},{pf_v},{pf_vw(end)},{pf_am(end)},{pf_Po(end)}];
smoother=PS_15(filter,a,o,C_noise);

ps=cell2mat(smoother(1));
ps_v=cell2mat(smoother(2));
ps_ml=transpose(ps(1,:));
ps_mr=transpose(ps(2,:));
ps_pl=transpose(ps(3,:));
ps_pr=transpose(ps(4,:));
ps_C =transpose(ps(5,:));

ps_rpl=1./(1+exp(-ps_ml));
ps_rpr=1./(1+exp(-ps_mr));
ps_precl=ps_pl./(ps_rpl.^2)./((1-ps_rpl).^2);
ps_precr=ps_pr./(ps_rpr.^2)./((1-ps_rpr).^2);

v_ml=ps_v(1,1,:);
V_ml=transpose(reshape(v_ml,1,length(a)));
v_mr=ps_v(2,2,:);
V_mr=transpose(reshape(v_mr,1,length(a)));
v_pl=ps_v(3,3,:);
V_pl=transpose(reshape(v_pl,1,length(a)));
v_pr=ps_v(4,4,:);
V_pr=transpose(reshape(v_pr,1,length(a)));
v_C=ps_v(5,5,:);
V_C=transpose(reshape(v_C,1,length(a)));

SD_rml=sqrt(V_ml.*(ps_rpl.^2).*(1-ps_rpl).^2);
SD_rmr=sqrt(V_mr.*(ps_rpr.^2).*(1-ps_rpr).^2);
SD_rpl=sqrt(V_pl./((ps_rpl.^4)./((1-ps_rpl).^4)));
SD_rpr=sqrt(V_pr./((ps_rpr.^4)./((1-ps_rpr).^4)));
SD_C=sqrt(V_C);

f = figure;
f.Position(3:4) = [1000 750];

N_step=length(a);

width=25;
for i=1:N_step-width
    ppl(i+width)=mean(a(i:i+width));
    
end

ppl(1:25)=NaN;
subplot(3,3,1:3);
plot(ppl)
xline(1,'-');xline(21,'-');xline(62,'-');
xline(124,'-');xline(218,'-');xline(238,'-');
xline(278,'-');xline(352,'-');
hold on;
 plot([0,length(a)],ones(1,2),'k')
     hold on
     plot([0,length(a)],0*ones(1,2),'k')
     ylabel('Right    Left','FontSize',15,'FontWeight','bold')
     yticks([0  0.5  1])
     ylim([-0.5 1.5])

 tmpr = find(a == 1 & o == 1);
 if ~isempty(tmpr), rtemp=plot(tmpr,1.25, 'b|','MarkerSize',30); 
     hold on; end
 
 tmp = find(a == 1 & o == 0);
 if ~isempty(tmp), ptemp=plot(tmp,1.25, 'b|','MarkerSize',15); 
     hold on; end
 tmpr = find(a == 0 & o == 1);
 if ~isempty(tmpr), rtemp=plot(tmpr,-0.25, 'r|','MarkerSize',30);
     hold on; end
 tmp = find(a == 0 & o == 0);
 if ~isempty(tmp), ptemp=plot(tmp,-0.25, 'r|','MarkerSize',15);
     hold on;end
     hold on;


subplot(3,3,4); 


plot(ppl)
hold on
plot(1-ppl)
xline(1,'-');xline(21,'-');xline(62,'-');
xline(124,'-');xline(218,'-');xline(238,'-');
xline(278,'-');xline(352,'-');
legend("left","right")
yticks([0 0.25 0.5 0.75 1])
ylim([0 1])


subplot(3,3,5); 
x=[1:1:length(a)];
sl1=0.9*ones(1,20);sl2=0.5*ones(1,41);sl3=0.5*ones(1,62);
sl4=0.1*ones(1,94);sl5=0.5*ones(1,20);sl6=0.9*ones(1,40);
sl7=0.1*ones(1,74);sl8=0.5*ones(1,40);
sl=[sl1,sl2,sl3,sl4,sl5,sl6,sl7,sl8];

sr1=0.5*ones(1,20);sr2=0.9*ones(1,41);sr3=0.1*ones(1,62);
sr4=0.5*ones(1,94);sr5=0.9*ones(1,20);sr6=0.5*ones(1,40);
sr7=0.5*ones(1,74);sr8=0.1*ones(1,40);
sr=[sr1,sr2,sr3,sr4,sr5,sr6,sr7,sr8];


shadedErrorBar(x,ps_rpl,SD_rml); 
hold on;
plot(x,ps_rpl,'- r','LineWidth',5)
hold on;
plot(x,sl,'-- r','LineWidth',3)
hold on;

xline(1,'-');xline(21,'-');xline(62,'-');
xline(124,'-');xline(218,'-');xline(238,'-');
xline(278,'-');xline(352,'-');
xlabel('Trials','FontSize',20,'FontWeight','bold')
yticks([0,0.25,0.5,0.75,1])
title('Reward Probability(Left)','FontSize',20,'FontWeight','bold')

subplot(3,3,6); 
x=[1:1:length(a)];
sl1=0.9*ones(1,20);sl2=0.5*ones(1,41);sl3=0.5*ones(1,62);
sl4=0.1*ones(1,94);sl5=0.5*ones(1,20);sl6=0.9*ones(1,40);
sl7=0.1*ones(1,74);sl8=0.5*ones(1,40);
sl=[sl1,sl2,sl3,sl4,sl5,sl6,sl7,sl8];

sr1=0.5*ones(1,20);sr2=0.9*ones(1,41);sr3=0.1*ones(1,62);
sr4=0.5*ones(1,94);sr5=0.9*ones(1,20);sr6=0.5*ones(1,40);
sr7=0.5*ones(1,74);sr8=0.1*ones(1,40);
sr=[sr1,sr2,sr3,sr4,sr5,sr6,sr7,sr8];


shadedErrorBar(x,ps_rpr,SD_rmr); 
hold on;
plot(x,ps_rpr,'- r','LineWidth',5)
hold on;
plot(x,sr,'-- r','LineWidth',3)
hold on;

xline(1,'-');xline(21,'-');xline(62,'-');
xline(124,'-');xline(218,'-');xline(238,'-');
xline(278,'-');xline(352,'-');

yticks([0,0.25,0.5,0.75,1])
xlabel('Trials','FontSize',20,'FontWeight','bold')
title('Reward Probability(Left)','FontSize',20,'FontWeight','bold')



subplot(3,3,7); 
shadedErrorBar(x,ps_precl,SD_rpl); 
hold on;
plot(ps_precl,'- r','LineWidth',5)
hold on;

xline(1,'-');xline(21,'-');xline(62,'-');
xline(124,'-');xline(218,'-');xline(238,'-');
xline(278,'-');xline(352,'-');
ylim([0 24])
yticks([0 6 12 18 24])
xlabel('Trials','FontSize',20,'FontWeight','bold')
title('Precision','FontSize',20,'FontWeight','bold')

subplot(3,3,8); 

shadedErrorBar(x,ps_precr,SD_rpr); 
hold on;
plot(ps_precr,'- g','LineWidth',5)
hold on;

xline(1,'-');xline(21,'-');xline(62,'-');
xline(124,'-');xline(218,'-');xline(238,'-');
xline(278,'-');xline(352,'-');
ylim([0 24])
yticks([0 6 12 18 24])
xlabel('Trials','FontSize',20,'FontWeight','bold')
title('Precision','FontSize',20,'FontWeight','bold')


subplot(3,3,9); 
hold on;
shadedErrorBar(x,ps_C,SD_C); 
plot(ps_C,'- r','LineWidth',5)
hold on;
xline(1,'-');xline(21,'-');xline(62,'-');
xline(124,'-');xline(218,'-');xline(238,'-');
xline(278,'-');xline(352,'-');
xlabel('Trials','FontSize',20,'FontWeight','bold')
title('Curiosity','FontSize',20','FontWeight','bold')
yline(0)
ylim([-9,3])
yticks([-9 -6 -3 0 3])
