%%
%Fig6 sample analysis code: effect of hyp/hss and targeted intervention
% Author: Chen Zhao, czhao22@jhmi.edu

%sensitivity analysis was done using the MATLAB code package provided in
%Marino et al. (PMID 18572196)


%load model
clear;
clc;
m = sbmlimport('7pathmodel_clean_v2.xml'); %import model: model name here should match the name of xml file
 cs = getconfigset(m, 'active');
cs.SolverOptions.AbsoluteTolerance=1e-9;
cs.SolverOptions.RelativeTolerance=1e-9;
cs.SolverType='ode15s';
 ics=[0	1200	468	17784	0	0	552129	78	28	36	58	36	1557	1	5425	10	1	4	1589	15960	7	87	1	2830	3	8	1	26	3068	0	0	59952	0	14	0	12	0	5	0	0	0	75355	18703	1	900	99122	93	1	16	0	1	1	1.45E+06	1	1	1	477	10	3	1	1	880	1437	1192	752	1	1	1	1227	873	14868	209	331	52	82	997	12376	1.204E+08	1405	35	8	384	30	5	5	1	1	1	1	151	1	1	1	1	17680	2	0	1	2	0	1	1	1	2450	6.50E+07	1	17649	249039	6	2	5	2	81	78	775	3828	173	131577	43860	18424	6141	1	4	759	0	0	6	0	1	4	24566	1	1	1	131	66	1	845	844	1	98868	1405	8	0	134379	886	1	3000	1	14869	132	675778	10093	14131	14521	144	100324	87	1	337	2000000	706	2470	5200	8	3	18	400	20	11	1	75892	4	1	277	931	25	29	68307	435	1	143	50	75239	4761	1	864	995	5	1	7100	1	995	5	1	845	1	625	1	199999	0	1	2995	146962	3120	5	11000	16	1	1	87270	52730];
newic=ics;

%%
%Compute relative fold change of M1-M2 markers under hyp/hss
%control
inos_b=newic(62);
arg1_b=newic(53);
socs3_b=newic(24)+newic(25)+newic(26)+newic(41)+newic(100);
c9_b=newic(83);
c10_b=newic(96);
il12_b=newic(90);
il1ra_b=newic(207);
mil10=newic(132);
nfkb=newic(176);
irf5=newic(205);
ap1=newic(187);
cebpb=newic(141);
h1=newic(75);
h2=newic(74);
ps3=newic(111);
ps1=newic(11);
ps6=newic(36);
cb=newic(185);
i1=newic(76);
i9=newic(107);
m93=newic(63);
il1bp_b=m.Parameters(261).Value*(nfkb*irf5*0.00033+ap1*cebpb*0.00001)/(nfkb*irf5*0.00033+ap1*cebpb*0.00001+30000)*(0.02+(h1*h2)^2/((h1*h2)^2+m.Parameters(284).Value))*(1.5-ps3/(ps3+m.Parameters(291).Value))*(1.05-(ps6^2)/(ps6^2+m.Parameters(300).Value));
ifngp_b=m.Parameters(30).Value*(0.01+(h1+h2)/(h1+h2+m.Parameters(296).Value))*(1-ps6/(ps6+m.Parameters(31).Value))*(0.1+(il12_b/(il12_b+m.Parameters(295).Value)))*(0.1+nfkb/(nfkb+m.Parameters(293).Value))*(0.1+(cb*ap1)/(cb*ap1+m.Parameters(294).Value))*(1-ps3/(ps3+m.Parameters(305).Value));
tnfap_b=m.Parameters(116).Value*(0.1+(i1^2)/(i1^2+m.Parameters(117).Value))*(15000+i9)*(nfkb*7*(irf5/(irf5+5000))^2+ap1*0.25+cebpb*0.0016-nfkb*cebpb*8.5e-8)/(nfkb*7*(irf5/(irf5+5000))^2+ap1*0.25+cebpb*0.0016-nfkb*cebpb*8.5e-8+100000)*(1.1-ps3/(ps3+m.Parameters(287).Value));
vap_b=m.Parameters(106).Value*(0.1+(h1*h2/(h1*h2+m.Parameters(107).Value)))*(1.1-m93/(m93+m.Parameters(310).Value))*(0.01+ap1/(ap1+m.Parameters(309).Value))*(cb/(cb+m.Parameters(316).Value));
il10p_b=m.Parameters(155).Value*mil10;


%Hyp/Hss
o2=1.204e8;
m.Species(78).InitialAmount=o2/21*2; % 2 percent o2
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 480);
[t1,out1,n2] = sbiosimulate(m);
inos=out1(end,62);
arg1=out1(end,53);
socs3=out1(end,24)+out1(end,25)+out1(end,26)+out1(end,41)+out1(end,100);
c9=out1(end,83);
c10=out1(end,96);
il12=out1(end,90);
il1ra=out1(end,207);
mil10=out1(end,132);
nfkb=out1(end,176);
irf5=out1(end,205);
ap1=out1(end,187);
cebpb=out1(end,141);
h1=out1(end,75);
h2=out1(end,74);
ps3=out1(end, 111);
ps1=out1(end,11);
ps6=out1(end,36);
cb=out1(end,185);
i1=out1(end,76);
i9=out1(end,107);
m93=out1(end,63);
il1bp=m.Parameters(261).Value*(nfkb*irf5*0.00033+ap1*cebpb*0.00001)/(nfkb*irf5*0.00033+ap1*cebpb*0.00001+30000)*(0.02+(h1*h2)^2/((h1*h2)^2+m.Parameters(284).Value))*(1.5-ps3/(ps3+m.Parameters(291).Value))*(1.05-(ps6^2)/(ps6^2+m.Parameters(300).Value));
ifngp=m.Parameters(30).Value*(0.01+(h1+h2)/(h1+h2+m.Parameters(296).Value))*(1-ps6/(ps6+m.Parameters(31).Value))*(0.1+(il12/(il12+m.Parameters(295).Value)))*(0.1+nfkb/(nfkb+m.Parameters(293).Value))*(0.1+(cb*ap1)/(cb*ap1+m.Parameters(294).Value))*(1-ps3/(ps3+m.Parameters(305).Value));
tnfap=m.Parameters(116).Value*(0.1+(i1^2)/(i1^2+m.Parameters(117).Value))*(15000+i9)*(nfkb*7*(irf5/(irf5+5000))^2+ap1*0.25+cebpb*0.0016-nfkb*cebpb*8.5e-8)/(nfkb*7*(irf5/(irf5+5000))^2+ap1*0.25+cebpb*0.0016-nfkb*cebpb*8.5e-8+100000)*(1.1-ps3/(ps3+m.Parameters(287).Value));
vap=m.Parameters(106).Value*(0.1+(h1*h2/(h1*h2+m.Parameters(107).Value)))*(1.1-m93/(m93+m.Parameters(310).Value))*(0.01+ap1/(ap1+m.Parameters(309).Value))*(cb/(cb+m.Parameters(316).Value));
il10p=m.Parameters(155).Value*mil10;

fc=[inos/inos_b arg1/arg1_b socs3/socs3_b c9/c9_b  tnfap/tnfap_b il1bp/il1bp_b il10p/il10p_b ifngp/ifngp_b c10/c10_b il12/il12_b vap/vap_b il1ra/il1ra_b ];
nfc=log2(fc);

%sample plot: Log2FC of M1-M2 markers at 8h of Hyp/hss, compare to Fig6B
figure(1);
bar([1 2 3 4 5 6 7 8 9 10 11 12],nfc);
%1-12: inos, arg1, socs3, cxcl9, tnfa, il1b, il10, ifng, cxcl10, il12,
%vegfa, il1ra
xtickangle(45);hold on;
axis([0 13 0 6]);
%%
%simulation of different targeted interventions under hyp/hss

m = sbmlimport('7pathmodel_clean_v2.xml'); %import model: model name here should match the name of xml file
 cs = getconfigset(m, 'active');
cs.SolverOptions.AbsoluteTolerance=1e-9;
cs.SolverOptions.RelativeTolerance=1e-9;
cs.SolverType='ode15s';

%overall M1/M2 score under control condition (2% O2)
m.Species(78).InitialAmount=1.204e8/21*2; %2 percent o2
set(cs, 'StopTime', 1500);
[t,out,n] = sbiosimulate(m);
inos=out(:,62)./ics(62);
arg1=out(:,53)./ics(53);
ifng=out(:,3)./ics(3);
tnfa=out(:,82)./ics(82);
il10=out(:,45)./ics(45);
il1b=out(:,196)./ics(196);
il1ra=out(:,207)./ics(207);
il12=out(:,90)./ics(90);
c9=out(:,83)./ics(83);
c10=out(:,96)./ics(96);
vegfa=out(:,2)./ics(2);
m1s=inos.*il12.*il1b.*tnfa.*ifng.*c9.*c10;
m2s=arg1.*il10.*vegfa.*il1ra;
m12=m1s./m2s;

%overall M1/M2 score under SOCS1 inhibition: SOCS1 degradation rate x10
m.Parameters(125).Value=0.004*10;
set(cs, 'StopTime', 1500);
[t2,out,n] = sbiosimulate(m);
inos=out(:,62)./ics(62);
arg1=out(:,53)./ics(53);
ifng=out(:,3)./ics(3);
tnfa=out(:,82)./ics(82);
il10=out(:,45)./ics(45);
il1b=out(:,196)./ics(196);
il1ra=out(:,207)./ics(207);
il12=out(:,90)./ics(90);
c9=out(:,83)./ics(83);
c10=out(:,96)./ics(96);
vegfa=out(:,2)./ics(2);
m1s_socs1=inos.*il12.*il1b.*tnfa.*ifng.*c9.*c10;
m2s_socs1=arg1.*il10.*vegfa.*il1ra;
m12_socs1=m1s_socs1./m2s_socs1;

%sample plot: change in relative M1, M2 and M1/M2 score under control and SOCS1 inhibition, compare to Fig6D
figure(2);
plot(t,log10(m1s),'m');
hold on;
plot(t2,log10(m1s_socs1),'m--');
plot(t,log10(m2s),'b');
hold on;
plot(t2,log10(m2s_socs1),'b--');
legend('M1-Ctr','M1-SOCS1 inhibition','M2-Ctr','M2-SOCS1 inhibition');
axis([0 1500 -10 10]);
figure(3);
plot(t,log10(m12),'r');
hold on;
plot(t2,log10(m12_socs1),'r--');
legend('M1/M2-Ctr','M1/M2-SOCS1 inhibition');
axis([0 1500 -10 10]);