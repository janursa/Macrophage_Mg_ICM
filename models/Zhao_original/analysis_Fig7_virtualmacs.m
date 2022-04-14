%%
%generation and analysis of virtual macrophage single cells
%Author: Chen Zhao, czhao22@jhmi.edu

%Note: Since the parameter sampling process (e.g. re-parameterization) is based on random number generator, for each run the set
%of 100 virtual single cells generated will give slightly different trajectories; the
%overall qualitative trends and conclusions from each run are consistent with the results presented in our paper. 

clear;
clc;
m = sbmlimport('7pathmodel_clean_v2.xml'); %import model: model name here should match the name of xml file
 cs = getconfigset(m, 'active');
cs.SolverOptions.AbsoluteTolerance=1e-9;
cs.SolverOptions.RelativeTolerance=1e-9;
cs.SolverType='ode15s';
 ics=[0	1200	468	17784	0	0	552129	78	28	36	58	36	1557	1	5425	10	1	4	1589	15960	7	87	1	2830	3	8	1	26	3068	0	0	59952	0	14	0	12	0	5	0	0	0	75355	18703	1	900	99122	93	1	16	0	1	1	1.45E+06	1	1	1	477	10	3	1	1	880	1437	1192	752	1	1	1	1227	873	14868	209	331	52	82	997	12376	1.204E+08	1405	35	8	384	30	5	5	1	1	1	1	151	1	1	1	1	17680	2	0	1	2	0	1	1	1	2450	6.50E+07	1	17649	249039	6	2	5	2	81	78	775	3828	173	131577	43860	18424	6141	1	4	759	0	0	6	0	1	4	24566	1	1	1	131	66	1	845	844	1	98868	1405	8	0	134379	886	1	3000	1	14869	132	675778	10093	14131	14521	144	100324	87	1	337	2000000	706	2470	5200	8	3	18	400	20	11	1	75892	4	1	277	931	25	29	68307	435	1	143	50	75239	4761	1	864	995	5	1	7100	1	995	5	1	845	1	625	1	199999	0	1	2995	146962	3120	5	11000	16	1	1	87270	52730];
%%
%resampling and re-parametrization to generate 100 new parameter sets that
%correspond to 100 virtual single cells
rng('shuffle');
plist=[254 187 70 69 94 95 93 272 277 278 269 318 319 320 321 137 117 287 116 284 291 300 261 296 31 295 293 294 305 30 288 154 289 290 302 303 118 107 310 309 316 106 ];
pps=zeros(100,42);
for i=1:100
p1=m.Parameters(254).Value; %ripk
r=2*rand-1;
k1=p1*(10^r);
k2=m.Parameters(187).Value*(3000/(3000+p1))/(3000/(3000+k1));
ps1dn=ics(11);
ps1i9=ics(22);
ps6dn=ics(36);

p3=m.Parameters(70).Value; %irf1
r=2*rand-1;
k3=p3*(10^r);
k4=m.Parameters(69).Value*(((ps1dn+ps1i9+1)/(ps6dn+1))/((ps1dn+ps1i9+1)/(ps6dn+1)+p3))/(((ps1dn+ps1i9+1)/(ps6dn+1))/((ps1dn+ps1i9+1)/(ps6dn+1)+p3));

pakt=ics(131); %irf4
p5=m.Parameters(94).Value;
r=2*rand-1;
k5=p5*(10^r);
p6=m.Parameters(95).Value;
r=2*rand-1;
k6=p6*(10^r);
k7=m.Parameters(93).Value*(ps6dn/(ps6dn+p5)*pakt/(pakt+p6))/(ps6dn/(ps6dn+k5)*pakt/(pakt+k6));

p8=m.Parameters(272).Value; %traf6
r=2*rand-1;
k8=p8*(10^r);
p9=m.Parameters(277).Value;
r=2*rand-1;
k9=p9*(10^r);
p10=m.Parameters(278).Value;
r=2*rand-1;
k10=p10*(10^r);
a20=ics(148);
s1=ics(15);
m146=ics(135);
k11=m.Parameters(269).Value*((1.01-a20/(a20+p8))*(1-s1/(s1+p9))*(1.2-m146/(m146+p10)))/((1.01-a20/(a20+k8))*(1-s1/(s1+k9))*(1.2-m146/(m146+k10)));

%socs1 
ps3dn=ics(111);
ap1=ics(187);
r=2*rand-1;
k12=83.5*(10^r);
r=2*rand-1;
k13=76*(10^r);
r=2*rand-1;
k14=2*(10^r);
r=2*rand-1;
k15=1000*(10^r);
k16=m.Parameters(137).Value*(0.1+(ps1dn/83.5)+(ps6dn/76)+(ps3dn/2)+(ap1/1000)^2)/((0.1+(ps1dn/k12)+(ps6dn/k13)+(ps3dn/k14)+(ap1/k15)^2));

%tnfa
i1=ics(76);
p17=m.Parameters(117).Value;
r=2*rand-1;
k17=p17*(10^r);
i9=ics(107);
nf=ics(176);
i5=ics(205);
ceb=ics(141);
p18=m.Parameters(287).Value;
r=2*rand-1;
k18=p18*(10^r);
k19=m.Parameters(116).Value*(0.1+(i1^2)/(i1^2+p17))*(15000+i9)*(nf*7*(i5/(i5+5000))^2+ap1*0.25+ceb*0.0016-nf*ceb*8.5e-8)/(nf*7*(i5/(i5+5000))^2+ap1*0.25+ceb*0.0016-nf*ceb*8.5e-8+100000)*(1.1-ps3dn/(ps3dn+p18))/((0.1+(i1^2)/(i1^2+k17))*(15000+i9)*(nf*7*(i5/(i5+5000))^2+ap1*0.25+ceb*0.0016-nf*ceb*8.5e-8)/(nf*7*(i5/(i5+5000))^2+ap1*0.25+ceb*0.0016-nf*ceb*8.5e-8+100000)*(1.1-ps3dn/(ps3dn+k18)));

%il1b
p20=m.Parameters(284).Value;
r=2*rand-1;
k20=p20*(10^r);
h1=ics(75);
h2=ics(74);
p21=m.Parameters(291).Value;
r=2*rand-1;
k21=p21*(10^r);
p22=m.Parameters(300).Value;
r=2*rand-1;
k22=p22*(10^r);
k23=m.Parameters(261).Value*(nf*i5*0.00033+ap1*ceb*0.00001)/(nf*i5*0.00033+ap1*ceb*0.00001+30000)*(0.02+(h1*h2)^2/((h1*h2)^2+p20))*(1.5-ps3dn/(ps3dn+p21))*(1.05-(ps6dn^2)/(ps6dn^2+p22))/((nf*i5*0.00033+ap1*ceb*0.00001)/(nf*i5*0.00033+ap1*ceb*0.00001+30000)*(0.02+(h1*h2)^2/((h1*h2)^2+k20))*(1.5-ps3dn/(ps3dn+k21))*(1.05-(ps6dn^2)/(ps6dn^2+k22)));

%IFNg
p24=m.Parameters(296).Value;
r=2*rand-1;
k24=p24*(10^r);
p25=m.Parameters(31).Value;
r=2*rand-1;
k25=p25*(10^r);
p26=m.Parameters(295).Value;
r=2*rand-1;
k26=p26*(10^r);
il12=ics(90);
p27=m.Parameters(293).Value;
r=2*rand-1;
k27=p27*(10^r);
cb=ics(187);
p28=m.Parameters(294).Value;
r=2*rand-1;
k28=p28*(10^r);
p29=m.Parameters(305).Value;
r=2*rand-1;
k29=p29*(10^r);
k30=m.Parameters(30).Value*(0.01+(h1+h2)/(h1+h2+p24))*(1-ps6dn/(ps6dn+p25))*(0.1+(il12/(il12+p26)))*(0.1+nf/(nf+p27))*(0.1+(cb*ap1)/(cb*ap1+p28))*(1-ps3dn/(ps3dn+p29))/((0.01+(h1+h2)/(h1+h2+k24))*(1-ps6dn/(ps6dn+k25))*(0.1+(il12/(il12+k26)))*(0.1+nf/(nf+k27))*(0.1+(cb*ap1)/(cb*ap1+k28))*(1-ps3dn/(ps3dn+k29)));

%il-10
p31=m.Parameters(288).Value;
r=2*rand-1;
k31=p31*(10^r);
p32=m.Parameters(154).Value;
r=2*rand-1;
k32=p32*(10^r);
p33=m.Parameters(289).Value;
r=2*rand-1;
k33=p33*(10^r);
p34=m.Parameters(290).Value;
r=2*rand-1;
k34=p34*(10^r);
p35=m.Parameters(302).Value;
r=2*rand-1;
k35=p35*(10^r);
p36=m.Parameters(303).Value;
r=2*rand-1;
k36=p36*(10^r);
k37=m.Parameters(118).Value*((ap1+cb*0.2+ceb*0.01)/(ap1+cb*0.2+ceb*0.01+p31))*(ps3dn)^1/((ps3dn)^1+p32)*(0.5+ps6dn/(ps6dn+p33))*(h1^2+p34)*(pakt/(pakt+p35))*(1-ps1dn/(ps1dn+p36))/(((ap1+cb*0.2+ceb*0.01)/(ap1+cb*0.2+ceb*0.01+k31))*(ps3dn)^1/((ps3dn)^1+k32)*(0.5+ps6dn/(ps6dn+k33))*(h1^2+k34)*(pakt/(pakt+k35))*(1-ps1dn/(ps1dn+k36)));

%vegfa
p38=m.Parameters(107).Value;
r=2*rand-1;
k38=p38*(10^r);
m93=ics(63);
p39=m.Parameters(310).Value;
r=2*rand-1;
k39=p39*(10^r);
p40=m.Parameters(309).Value;
r=2*rand-1;
k40=p40*(10^r);
p41=m.Parameters(316).Value;
r=2*rand-1;
k41=p41*(10^r);
k42=m.Parameters(106).Value*(0.1+(h1*h2/(h1*h2+p38)))*(1.1-m93/(m93+p39))*(0.01+ap1/(ap1+p40))*(cb/(cb+p41))/((0.1+(h1*h2/(h1*h2+k38)))*(1.1-m93/(m93+k39))*(0.01+ap1/(ap1+k40))*(cb/(cb+k41)));

pps(i,1:42)=[k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17 k18 k19 k20 k21 k22 k23 k24 k25 k26 k27 k28 k29 k30 k31 k32 k33 k34 k35 k36 k37 k38 k39 k40 k41 k42];
end
%%
%simulate and plot the distribution of M1, M2, M1/M2 response (in 100 single cells) under
%control condition (2% O2)
m = sbmlimport('7pathmodel_clean_v2.xml'); %import model: model name here should match the name of xml file
 cs = getconfigset(m, 'active');
cs.SolverOptions.AbsoluteTolerance=1e-9;
cs.SolverOptions.RelativeTolerance=1e-9;
cs.SolverType='ode15s';
m.Species(78).InitialAmount=1.204e8/21*2; 
m1c=0;
m2c=0;
m0c=0;
for oi=1:100;
for i=1:42
    m.Parameters(plist(i)).Value=pps(oi,i);
end

set(cs, 'StopTime', 1500);
[t,out] = sbiosimulate(m);
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
 if (log10(m12(end))>2)
    m1c=m1c+1;
    cntm1(m1c)=oi;
  elseif (log10(m12(end))<-1)
    m2c=m2c+1;
    cntm2(m2c)=oi;
 else 
    m0c=m0c+1;
    cntm0(m0c)=oi;
end
figure(1);
aa=plot(t,log10(m1s),'LineWidth',2);
aa.Color(4)=0.4; hold on;
figure(2);
aa=plot(t,log10(m2s),'LineWidth',2);
aa.Color(4)=0.4; hold on;
figure(3);
aa=plot(t,log10(m12),'LineWidth',2);
aa.Color(4)=0.4; hold on;
end

m = sbmlimport('7pathmodel_clean_v2.xml'); %import model: model name here should match the name of xml file
 cs = getconfigset(m, 'active');
cs.SolverOptions.AbsoluteTolerance=1e-9;
cs.SolverOptions.RelativeTolerance=1e-9;
cs.SolverType='ode15s';
m.Species(78).InitialAmount=1.204e8/21*2; 

set(cs, 'StopTime', 1500);
[t,out] = sbiosimulate(m);
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


%sample plots: compare to Figs7A-C
figure(1);
plot(t, log10(m1s),'k','LineWidth',4);hold on;
axis([0 1500 -5 10]);
figure(2);
plot(t, log10(m2s),'k','LineWidth',4);hold on;
axis([0 1500 -5 10]);
figure(3);
plot(t, log10(m12),'k','LineWidth',4);hold on;
axis([0 1500 -10 10]);
%%
%simulate and plot the distribution of M1/M2 response (in 100 single cells) under
%control condition (2% O2) with an additional perturbation
m = sbmlimport('7pathmodel_clean_v2.xml'); %import model: model name here should match the name of xml file
 cs = getconfigset(m, 'active');
cs.SolverOptions.AbsoluteTolerance=1e-9;
cs.SolverOptions.RelativeTolerance=1e-9;
cs.SolverType='ode15s';
m.Species(78).InitialAmount=1.204e8/21*2; 
m1c=0;
m2c=0;
m0c=0;
figure(4);
for oi=1:100;
for i=1:42
    m.Parameters(plist(i)).Value=pps(oi,i);
end
m.Parameters(28).Value=0.01*0.1; %pstat6 deact rate x0.1
%m.Parameters(129).Value=0.3*0.1; %pstat3 deact rate x0.1
set(cs, 'StopTime', 1500);
[t,out] = sbiosimulate(m);
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
 if (log10(m12(end))>2)
    m1c=m1c+1;
    cntm1(m1c)=oi;
  elseif (log10(m12(end))<-1)
    m2c=m2c+1;
    cntm2(m2c)=oi;
 else 
    m0c=m0c+1;
    cntm0(m0c)=oi;
end

aa=plot(t,log10(m12),'LineWidth',2);
aa.Color(4)=0.4; hold on;
end

m = sbmlimport('7pathmodel_clean_v2.xml'); %import model: model name here should match the name of xml file
 cs = getconfigset(m, 'active');
cs.SolverOptions.AbsoluteTolerance=1e-9;
cs.SolverOptions.RelativeTolerance=1e-9;
cs.SolverType='ode15s';
m.Species(78).InitialAmount=1.204e8/21*2; 
m.Parameters(28).Value=0.01*0.1; %pstat6 deact rate x0.1
%m.Parameters(129).Value=0.3*0.1; %pstat3 deact rate x0.1
set(cs, 'StopTime', 1500);
[t,out] = sbiosimulate(m);
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

%sample plot: compare to Fig7D (STAT6*) and Fig7E (STAT3*)
figure(4);
plot(t, log10(m12),'k','LineWidth',4);hold on;
axis([0 1500 -15 5]); %stat6* axis range
%axis([0 1500 -10 10]); %stat3* axis range