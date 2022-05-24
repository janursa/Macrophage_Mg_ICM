clear;
clc;
m = sbmlimport('7pathmodel_clean_v2.xml'); %import model: model name here should match the name of xml file
species = m.Species;
cs = getconfigset(m, 'active');
 
cs.SolverOptions.AbsoluteTolerance=1e-9;
cs.SolverOptions.RelativeTolerance=1e-9;
cs.SolverType='ode15s';
set(cs, 'StopTime', 1500); %set simulation time (in minutes)
[t,out] = sbiosimulate(m); %simulate model
h1 = out(:,82);
figure(1); %reproduce the simulated profile of HIF1a in Fig.2Y
plot(t,h1,'r');
hold on;
xlabel('Time (min)');
ylabel('TNFa');