%%
%Fig5 sample analysis code: generating values for heatmap and time-course
%M1/M2 score curves (in vitro scenario)
% Author: Chen Zhao, czhao22@jhmi.edu


clear;
clc;
m = sbmlimport('7pathmodel_clean_v2.xml'); %import model: model name here should match the name of xml file
 cs = getconfigset(m, 'active');
cs.SolverOptions.AbsoluteTolerance=1e-9;
cs.SolverOptions.RelativeTolerance=1e-9;
cs.SolverType='ode15s';
 ics=[0	1200	468	17784	0	0	552129	78	28	36	58	36	1557	1	5425	10	1	4	1589	15960	7	87	1	2830	3	8	1	26	3068	0	0	59952	0	14	0	12	0	5	0	0	0	75355	18703	1	900	99122	93	1	16	0	1	1	1.45E+06	1	1	1	477	10	3	1	1	880	1437	1192	752	1	1	1	1227	873	14868	209	331	52	82	997	12376	1.204E+08	1405	35	8	384	30	5	5	1	1	1	1	151	1	1	1	1	17680	2	0	1	2	0	1	1	1	2450	6.50E+07	1	17649	249039	6	2	5	2	81	78	775	3828	173	131577	43860	18424	6141	1	4	759	0	0	6	0	1	4	24566	1	1	1	131	66	1	845	844	1	98868	1405	8	0	134379	886	1	3000	1	14869	132	675778	10093	14131	14521	144	100324	87	1	337	2000000	706	2470	5200	8	3	18	400	20	11	1	75892	4	1	277	931	25	29	68307	435	1	143	50	75239	4761	1	864	995	5	1	7100	1	995	5	1	845	1	625	1	199999	0	1	2995	146962	3120	5	11000	16	1	1	87270	52730];
newic=ics;

%values in control condition
inos_b=newic(62);
arg_b=newic(53);
c9_b=newic(83);
c10_b=newic(96);
il12_b=newic(90);
il1ra_b=newic(207);
mil10_b=newic(132);
nfkb_b=newic(176);
irf5_b=newic(205);
ap1_b=newic(187);
cebpb_b=newic(141);
h1_b=newic(75);
h2_b=newic(74);
ps3_b=newic(111);
ps1_b=newic(11);
ps6_b=newic(36);
cb_b=newic(185);
i1_b=newic(76);
i9_b=newic(107);
m93_b=newic(63);
il1bp_b=m.Parameters(261).Value*(nfkb_b*irf5_b*0.00033+ap1_b*cebpb_b*0.00001)/(nfkb_b*irf5_b*0.00033+ap1_b*cebpb_b*0.00001+30000)*(0.02+(h1_b*h2_b)^2/((h1_b*h2_b)^2+m.Parameters(284).Value))*(1.5-ps3_b/(ps3_b+m.Parameters(291).Value))*(1.05-(ps6_b^2)/(ps6_b^2+m.Parameters(300).Value));
ifngp_b=m.Parameters(30).Value*(0.01+(h1_b+h2_b)/(h1_b+h2_b+m.Parameters(296).Value))*(1-ps6_b/(ps6_b+m.Parameters(31).Value))*(0.1+(il12_b/(il12_b+m.Parameters(295).Value)))*(0.1+nfkb_b/(nfkb_b+m.Parameters(293).Value))*(0.1+(cb_b*ap1_b)/(cb_b*ap1_b+m.Parameters(294).Value))*(1-ps3_b/(ps3_b+m.Parameters(305).Value));
tnfap_b=m.Parameters(116).Value*(0.1+(i1_b^2)/(i1_b^2+m.Parameters(117).Value))*(15000+i9_b)*(nfkb_b*7*(irf5_b/(irf5_b+5000))^2+ap1_b*0.25+cebpb_b*0.0016-nfkb_b*cebpb_b*8.5e-8)/(nfkb_b*7*(irf5_b/(irf5_b+5000))^2+ap1_b*0.25+cebpb_b*0.0016-nfkb_b*cebpb_b*8.5e-8+100000)*(1.1-ps3_b/(ps3_b+m.Parameters(287).Value));
vap_b=m.Parameters(106).Value*(0.1+(h1_b*h2_b/(h1_b*h2_b+m.Parameters(107).Value)))*(1.1-m93_b/(m93_b+m.Parameters(310).Value))*(0.01+ap1_b/(ap1_b+m.Parameters(309).Value))*(cb_b/(cb_b+m.Parameters(316).Value));
il10p_b=m.Parameters(155).Value*mil10_b;
base=[inos_b arg_b c9_b c10_b il12_b il1ra_b il1bp_b ifngp_b tnfap_b vap_b il10p_b];
sb=inos_b*c9_b*c10_b*il12_b*il1bp_b*ifngp_b*tnfap_b/arg_b/il1ra_b/vap_b/il10p_b;

%simulating treatment combinations (assuming in vitro doses): at 4h, 24h, 48h
%hypoxia is used as an example

%hyp 2% o2 4h
o2=1.204e8;
m.Species(78).InitialAmount=o2/21*2; % 2 percent o2
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hyp_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hyp_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;



%hyp 2% o2 24h

m.Species(78).InitialAmount=o2/21*2; % 2 percent o2
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hyp_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hyp_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp 2% o2 48h

m.Species(78).InitialAmount=o2/21*2; % 2 percent o2
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hyp_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hyp_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
hyp_s=[1 hyp_4s/sb hyp_24s/sb hyp_48s/sb];

%tnfa 10 ng/ml: 4h

m.Species(78).InitialAmount=newic(78);
m.Species(82).InitialAmount=110000; %10 ng/ml tNFa
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnf_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnf_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%tnfa 10 ng/ml: 24h

m.Species(78).InitialAmount=newic(78);
m.Species(82).InitialAmount=110000; %10 ng/ml tNFa
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnf_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnf_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%tnfa 10 ng/ml: 48h

m.Species(78).InitialAmount=newic(78);
m.Species(82).InitialAmount=110000; %10 ng/ml tNFa
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnf_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnf_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
tnf_s=[1 tnf_4s/sb tnf_24s/sb tnf_48s/sb];

%ifng 10 ng/ml: 4h

m.Species(82).InitialAmount=newic(82);
m.Species(3).InitialAmount=300000; %10 ng/ml IFNg
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifn_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifn_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng 10 ng/ml: 24h

m.Species(82).InitialAmount=newic(82);
m.Species(3).InitialAmount=300000; %10 ng/ml IFNg
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifn_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifn_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng 10 ng/ml: 48h

m.Species(82).InitialAmount=newic(82);
m.Species(3).InitialAmount=300000; %10 ng/ml IFNg
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;
ifn_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifn_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
ifn_s=[1 ifn_4s/sb ifn_24s/sb ifn_48s/sb];

%il1b 10 ng/ml: 4h

m.Species(3).InitialAmount=newic(3);
m.Species(196).InitialAmount=350000; %10 ng/ml IL1b
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;
il1_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il1b 10 ng/ml: 24h

m.Species(3).InitialAmount=newic(3);
m.Species(196).InitialAmount=350000; %10 ng/ml IL1b
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il1_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il1b 10 ng/ml: 48h

m.Species(3).InitialAmount=newic(3);
m.Species(196).InitialAmount=350000; %10 ng/ml IL1b
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;
il1_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
il1_s=[1 il1_4s/sb il1_24s/sb il1_48s/sb];

%il-4 10 ng/ml: 4h

m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000; %10 ng/ml IL-4
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;
il4_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il4_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il-4 10 ng/ml: 24h

m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000; %10 ng/ml IL-4
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;
il4_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il4_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il-4 10 ng/ml: 48h

m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000; %10 ng/ml IL-4
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il4_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il4_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
il4_s=[1 il4_4s/sb il4_24s/sb il4_48s/sb];

%il-10 10 ng/ml: 4h

m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000; %10 ng/ml IL-10
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il10_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il10_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il-10 10 ng/ml: 24h

m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000; %10 ng/ml IL-10
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);
il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il10_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il10_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il-10 10 ng/ml: 48h

m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000; %10 ng/ml IL-10
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il10_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il10_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
il10_s=[1 il10_4s/sb il10_24s/sb il10_48s/sb];

%vegf 10 ng/ml: 4h

m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000; %10 ng/ml VEGF
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

vegf_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
vegf_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%vegf 10 ng/ml: 24h

m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000; %10 ng/ml VEGF
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

vegf_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
vegf_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%vegf 10 ng/ml: 48h

m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000; %10 ng/ml VEGF
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

vegf_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
vegf_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
vegf_s=[1 vegf_4s/sb vegf_24s/sb vegf_48s/sb];

%hyp+ifn 4h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=300000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypifn_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypifn_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+ifn 24 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=300000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypifn_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypifn_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+ifn 48 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=300000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypifn_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypifn_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
hypifn_s=[1 hypifn_4s/sb hypifn_24s/sb hypifn_48s/sb];

%hyp+tnf 4h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hyptnf_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hyptnf_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+tnf 24 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hyptnf_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hyptnf_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+tnf 48 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;
hyptnf_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hyptnf_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
hyptnf_s=[1 hyptnf_4s/sb hyptnf_24s/sb hyptnf_48s/sb];

%hyp+il1 4 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000; %10 ng/ml IL-1
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypil1_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypil1_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+il1 24h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000; %10 ng/ml IL-1
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypil1_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypil1_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+il1 48 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000; %10 ng/ml IL-1
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypil1_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypil1_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
hypil1_s=[1 hypil1_4s/sb hypil1_24s/sb hypil1_48s/sb];

%hyp+IL-4 4 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000; %10 ng/ml il-4
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypil4_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypil4_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+IL-4 24 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000; %10 ng/ml il-4
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypil4_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypil4_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+IL-4 48 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000; %10 ng/ml il-4
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypil4_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypil4_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
hypil4_s=[1 hypil4_4s/sb hypil4_24s/sb hypil4_48s/sb];

%hyp+IL-10 4h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypil10_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypil10_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+IL-10 24h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypil10_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypil10_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+IL-10 48 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypil10_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypil10_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
hypil10_s=[1 hypil10_4s/sb hypil10_24s/sb hypil10_48s/sb];

%hyp+vegf 4 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000; %10 ng/ml VEGF
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;
hypvegf_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypvegf_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+vegf 24 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000; %10 ng/ml VEGF
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypvegf_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypvegf_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%hyp+vegf 48 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=o2/21*2;
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000; %10 ng/ml VEGF
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

hypvegf_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
hypvegf_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
hypvegf_s=[1 hypvegf_4s/sb hypvegf_24s/sb hypvegf_48s/sb];

%ifng+tnfa 4 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifntnf_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifntnf_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng+tnfa 24 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifntnf_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifntnf_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng+tnfa 48 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifntnf_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifntnf_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
ifntnf_s=[1 ifntnf_4s/sb ifntnf_24s/sb ifntnf_48s/sb];

%ifng+il1 4 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnil1_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnil1_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng+il1 24 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnil1_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnil1_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng+il1 48 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnil1_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnil1_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
ifnil1_s=[1 ifnil1_4s/sb ifnil1_24s/sb ifnil1_48s/sb];

%ifng+il4 4 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnil4_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnil4_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng+il4 24 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnil4_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnil4_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng+il4 48 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnil4_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnil4_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
ifnil4_s=[1 ifnil4_4s/sb ifnil4_24s/sb ifnil4_48s/sb];

%ifng+il10 4 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnil10_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnil10_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng+il10 24 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnil10_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnil10_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng+il10 48 h

m.Species(2).InitialAmount=newic(2);
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnil10_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnil10_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
ifnil10_s=[1 ifnil10_4s/sb ifnil10_24s/sb ifnil10_48s/sb];

%ifng+vegf 4h

m.Species(2).InitialAmount=290000;
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnvegf_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnvegf_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng+vegf 24 h

m.Species(2).InitialAmount=290000;
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnvegf_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnvegf_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%ifng+vegf 48 h

m.Species(2).InitialAmount=290000;
m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=300000;
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

ifnvegf_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
ifnvegf_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
ifnvegf_s=[1 ifnvegf_4s/sb ifnvegf_24s/sb ifnvegf_48s/sb];

%tnfa+il1 4 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfil1_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfil1_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%tnfa+il1 24 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfil1_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfil1_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%tnfa+il1 48 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfil1_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfil1_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
tnfil1_s=[1 tnfil1_4s/sb tnfil1_24s/sb tnfil1_48s/sb];

%tnfa+il4 4 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfil4_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfil4_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%tnfa+il4 24 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfil4_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfil4_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%tnfa+il4 48 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfil4_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfil4_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
tnfil4_s=[1 tnfil4_4s/sb tnfil4_24s/sb tnfil4_48s/sb];

%tnfa+il10 4 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfil10_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfil10_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%tnfa+il10 24 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfil10_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfil10_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%tnfa+il10 48 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfil10_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfil10_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
tnfil10_s=[1 tnfil10_4s/sb tnfil10_24s/sb tnfil10_48s/sb];

%tnfa+vegf 4h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfvegf_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfvegf_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%tnfa+vegf 24 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfvegf_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfvegf_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%tnfa+vegf 48 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=110000;
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

tnfvegf_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
tnfvegf_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
tnfvegf_s=[1 tnfvegf_4s/sb tnfvegf_24s/sb tnfvegf_48s/sb];

%il1+il4 4h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il1il4_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1il4_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il1+il4 24 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il1il4_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1il4_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il1+il4 48 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il1il4_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1il4_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
il1il4_s=[1 il1il4_4s/sb il1il4_24s/sb il1il4_48s/sb];

%il1+il10 4h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il1il10_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1il10_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il1+il10 24 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il1il10_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1il10_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il1+il10 48 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il1il10_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1il10_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
il1il10_s=[1 il1il10_4s/sb il1il10_24s/sb il1il10_48s/sb];

%il1+vegf 4 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il1vegf_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1vegf_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il1+vegf 24 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il1vegf_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1vegf_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il1+vegf 48 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=350000;
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il1vegf_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il1vegf_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
il1vegf_s=[1 il1vegf_4s/sb il1vegf_24s/sb il1vegf_48s/sb];

%il4+il10 4 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il4il10_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il4il10_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il4+il10 24 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il4il10_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il4il10_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il4+il10 48 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=newic(2);
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il4il10_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il4il10_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
il4il10_s=[1 il4il10_4s/sb il4il10_24s/sb il4il10_48s/sb];

%il4+vegf 4 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il4vegf_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il4vegf_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il4+vegf 24 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il4vegf_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il4vegf_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il4+vegf 48 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=400000;
m.Species(45).InitialAmount=newic(45);
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il4vegf_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il4vegf_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
il4vegf_s=[1 il4vegf_4s/sb il4vegf_24s/sb il4vegf_48s/sb];

%il10+vegf 4 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 240);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il10vegf_4=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il10vegf_4s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il10+vegf 24 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 1440);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;
il10vegf_24=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il10vegf_24s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;

%il10+vegf 48 h

m.Species(78).InitialAmount=newic(78);
m.Species(3).InitialAmount=newic(3);
m.Species(82).InitialAmount=newic(82);
m.Species(196).InitialAmount=newic(196);
m.Species(28).InitialAmount=newic(28);
m.Species(45).InitialAmount=300000;
m.Species(2).InitialAmount=290000;
cs = getconfigset(m, 'active');
set(cs, 'StopTime', 2880);
[t1,out1,n2] = sbiosimulate(m);
inos_sti=out1(end,62);
arg_sti=out1(end,53);
c9_sti=out1(end,83);
c10_sti=out1(end,96);
il12_sti=out1(end,90);
il1ra_sti=out1(end,207);
mil10_sti=out1(end,132);
nfkb_sti=out1(end,176);
irf5_sti=out1(end,205);
ap1_sti=out1(end,187);
cebpb_sti=out1(end,141);
h1_sti=out1(end,75);
h2_sti=out1(end,74);
ps3_sti=out1(end, 111);
ps1_sti=out1(end,11);
ps6_sti=out1(end,36);
cb_sti=out1(end,185);
i1_sti=out1(end,76);
i9_sti=out1(end,107);
m93_sti=out1(end,63);

il1bp_sti=m.Parameters(261).Value*(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001)/(nfkb_sti*irf5_sti*0.00033+ap1_sti*cebpb_sti*0.00001+30000)*(0.02+(h1_sti*h2_sti)^2/((h1_sti*h2_sti)^2+m.Parameters(284).Value))*(1.5-ps3_sti/(ps3_sti+m.Parameters(291).Value))*(1.05-(ps6_sti^2)/(ps6_sti^2+m.Parameters(300).Value));
ifngp_sti=m.Parameters(30).Value*(0.01+(h1_sti+h2_sti)/(h1_sti+h2_sti+m.Parameters(296).Value))*(1-ps6_sti/(ps6_sti+m.Parameters(31).Value))*(0.1+(il12_sti/(il12_sti+m.Parameters(295).Value)))*(0.1+nfkb_sti/(nfkb_sti+m.Parameters(293).Value))*(0.1+(cb_sti*ap1_sti)/(cb_sti*ap1_sti+m.Parameters(294).Value))*(1-ps3_sti/(ps3_sti+m.Parameters(305).Value));
tnfap_sti=m.Parameters(116).Value*(0.1+(i1_sti^2)/(i1_sti^2+m.Parameters(117).Value))*(15000+i9_sti)*(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8)/(nfkb_sti*7*(irf5_sti/(irf5_sti+5000))^2+ap1_sti*0.25+cebpb_sti*0.0016-nfkb_sti*cebpb_sti*8.5e-8+100000)*(1.1-ps3_sti/(ps3_sti+m.Parameters(287).Value));
vap_sti=m.Parameters(106).Value*(0.1+(h1_sti*h2_sti/(h1_sti*h2_sti+m.Parameters(107).Value)))*(1.1-m93_sti/(m93_sti+m.Parameters(310).Value))*(0.01+ap1_sti/(ap1_sti+m.Parameters(309).Value))*(cb_sti/(cb_sti+m.Parameters(316).Value));
il10p_sti=m.Parameters(155).Value*mil10_sti;

il10vegf_48=[inos_sti arg_sti c9_sti c10_sti il12_sti il1ra_sti il1bp_sti ifngp_sti tnfap_sti vap_sti il10p_sti];
il10vegf_48s=inos_sti*c9_sti*c10_sti*il12_sti*il1bp_sti*ifngp_sti*tnfap_sti/arg_sti/il1ra_sti/vap_sti/il10p_sti;
il10vegf_s=[1 il10vegf_4s/sb il10vegf_24s/sb il10vegf_48s/sb];

%4,24 ,48h
%these are data value inputs for the heatmap as shown in Fig5A
sum4vitro=[hyp_4./base;hypifn_4./base; hyptnf_4./base; hypil1_4./base; hypil4_4./base; hypil10_4./base; hypvegf_4./base; ifn_4./base; ifntnf_4./base; ifnil1_4./base; ifnil4_4./base; ifnil10_4./base; ifnvegf_4./base; tnf_4./base; tnfil1_4./base; tnfil4_4./base; tnfil10_4./base; tnfvegf_4./base;il1_4./base; il1il4_4./base; il1il10_4./base; il1vegf_4./base; il4_4./base; il4il10_4./base; il4vegf_4./base; il10_4./base; il10vegf_4./base; vegf_4./base];
sum24vitro=[hyp_24./base;hypifn_24./base; hyptnf_24./base; hypil1_24./base; hypil4_24./base; hypil10_24./base; hypvegf_24./base; ifn_24./base; ifntnf_24./base; ifnil1_24./base; ifnil4_24./base; ifnil10_24./base; ifnvegf_24./base; tnf_24./base; tnfil1_24./base; tnfil4_24./base; tnfil10_24./base; tnfvegf_24./base;il1_24./base; il1il4_24./base; il1il10_24./base; il1vegf_24./base; il4_24./base; il4il10_24./base; il4vegf_24./base; il10_24./base; il10vegf_24./base; vegf_24./base];
sum48vitro=[hyp_48./base;hypifn_48./base; hyptnf_48./base; hypil1_48./base; hypil4_48./base; hypil10_48./base; hypvegf_48./base; ifn_48./base; ifntnf_48./base; ifnil1_48./base; ifnil4_48./base; ifnil10_48./base; ifnvegf_48./base; tnf_48./base; tnfil1_48./base; tnfil4_48./base; tnfil10_48./base; tnfvegf_48./base;il1_48./base; il1il4_48./base; il1il10_48./base; il1vegf_48./base; il4_48./base; il4il10_48./base; il4vegf_48./base; il10_48./base; il10vegf_48./base; vegf_48./base];
%%
%sample plot: time course M1/M2 score trends for all conditions (in vitro)
%compare to Fig 5B (8 curves are bolded, with different color coding)
figure(1); %in vitro trend changes
plot([0 240 1440 2880],hyp_s);
hold on;
plot([0 240 1440 2880],hypifn_s,'LineWidth',3);
plot([0 240 1440 2880],hyptnf_s);
plot([0 240 1440 2880],hypil1_s);
plot([0 240 1440 2880],hypil4_s);
plot([0 240 1440 2880],hypil10_s);
plot([0 240 1440 2880],hypvegf_s);
plot([0 240 1440 2880],ifn_s,'LineWidth',3);
plot([0 240 1440 2880],ifntnf_s,'LineWidth',3);
plot([0 240 1440 2880],ifnil1_s);
plot([0 240 1440 2880],ifnil4_s,'LineWidth',3);
plot([0 240 1440 2880],ifnil10_s);
plot([0 240 1440 2880],ifnvegf_s);
plot([0 240 1440 2880],tnf_s);
plot([0 240 1440 2880],tnfil1_s);
plot([0 240 1440 2880],tnfil4_s);
plot([0 240 1440 2880],tnfil10_s);
plot([0 240 1440 2880],tnfvegf_s);
plot([0 240 1440 2880],il1_s);
plot([0 240 1440 2880],il1il4_s);
plot([0 240 1440 2880],il1il10_s,'LineWidth',3);
plot([0 240 1440 2880],il1vegf_s);
plot([0 240 1440 2880],il4_s,'LineWidth',3);
plot([0 240 1440 2880],il4il10_s,'LineWidth',3);
plot([0 240 1440 2880],il4vegf_s);
plot([0 240 1440 2880],il10_s,'LineWidth',3);
plot([0 240 1440 2880],il10vegf_s);
plot([0 240 1440 2880],vegf_s);
set(gca, 'YScale', 'log')
axis([0 2880 1e-13 1e11]);