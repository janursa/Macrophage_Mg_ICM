%%
%model setup and simulation
% Author: Chen Zhao, czhao22@jhmi.edu

clear;
clc;
m = sbmlimport('edited.xml'); %import model: model name here should match the name of xml file
 cs = getconfigset(m, 'active');
cs.SolverOptions.AbsoluteTolerance=1e-9;
cs.SolverOptions.RelativeTolerance=1e-9;
cs.SolverType='ode15s';
 ics=[0	1200	468	17784	0	0	552129	78	28	36	58	36	1557	1	5425	10	1	4	1589	15960	7	87	1	2830	3	8	1	26	3068	0	0	59952	0	14	0	12	0	5	0	0	0	75355	18703	1	900	99122	93	1	16	0	1	1	1.45E+06	1	1	1	477	10	3	1	1	880	1437	1192	752	1	1	1	1227	873	14868	209	331	52	82	997	12376	1.204E+08	1405	35	8	384	30	5	5	1	1	1	1	151	1	1	1	1	17680	2	0	1	2	0	1	1	1	2450	6.50E+07	1	17649	249039	6	2	5	2	81	78	775	3828	173	131577	43860	18424	6141	1	4	759	0	0	6	0	1	4	24566	1	1	1	131	66	1	845	844	1	98868	1405	8	0	134379	886	1	3000	1	14869	132	675778	10093	14131	14521	144	100324	87	1	337	2000000	706	2470	5200	8	3	18	400	20	11	1	75892	4	1	277	931	25	29	68307	435	1	143	50	75239	4761	1	864	995	5	1	7100	1	995	5	1	845	1	625	1	199999	0	1	2995	146962	3120	5	11000	16	1	1	87270	52730];
species = m.Species();



%% 
%simulation scenario 1: hypoxia
[t,out]=sbiosimulate(m);
o2=ics(78);
m.Species(78).InitialAmount=o2/21*3; %simulate 3 percent oxyge
set(cs, 'StopTime', 1500); %set simulation time (in minutes)
[t,out] = sbiosimulate(m); %simulate model


h1=out(:,64)+out(:,73)+out(:,75)+out(:,80); %get simulated total cellular hif1a expression
figure(1); %reproduce the simulated profile of HIF1a in Fig.2Y
plot(t,h1./max(h1),'r');
hold on;
xlabel('Time (min)');
ylabel('Relative Expression');
legend('Hypoxia induces HIF1a expression');
axis([0 1500 0 1.2]);

%% 
%simulation scenario 2: IL-4
m.Species(78).InitialAmount=ics(78);%restore IC for oxygen
m.Species(28)
m.Species(28).InitialAmount=400000; %simulate 10 ng/ml IL-4
cs = getconfigset(m, 'active');
 set(cs, 'StopTime', 1500); 
[t,out] = sbiosimulate(m);
a1=out(:,45); %get simulated arg1 expression 
figure(2); %reproduce the simulated profile of ARG1 in Fig.3G
plot(t,a1./max(a1),'r');
hold on;
xlabel('Time (min)');
ylabel('Relative Expression');
legend('IL-10 ');
axis([0 1500 0 1.2]);
