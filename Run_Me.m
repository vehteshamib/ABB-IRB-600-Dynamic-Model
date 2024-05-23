clear all
close all
clc

global g t1 t2 t3 t4 t5 t6   
global xJ1 yJ1 zJ1   
global xG1 yG1 zG1 m1 Ixx1 Iyy1 Izz1 Ixy1 Ixz1 Iyz1   
global xJ2 yJ2 zJ2   
global xG2 yG2 zG2 m2 Ixx2 Iyy2 Izz2 Ixy2 Ixz2 Iyz2   
global xJ3 yJ3 zJ3   
global xG3 yG3 zG3 m3 Ixx3 Iyy3 Izz3 Ixy3 Ixz3 Iyz3   
global xJ4 yJ4 zJ4   
global xG4 yG4 zG4 m4 Ixx4 Iyy4 Izz4 Ixy4 Ixz4 Iyz4   
global xJ5 yJ5 zJ5   
global xG5 yG5 zG5 m5 Ixx5 Iyy5 Izz5 Ixy5 Ixz5 Iyz5   
global xJ6 yJ6 zJ6   
global xG6 yG6 zG6 m6 Ixx6 Iyy6 Izz6 Ixy6 Ixz6 Iyz6 

% Geometric Parameters
 xJ1=0; yJ1=0; zJ1=0;
 xG1=0; yG1=0.100; zG1=0.250;
 m1=20; Ixx1=0.2; Iyy1=0.2; Izz1=0.2; Ixy1=0; Ixz1=0; Iyz1=0;   
 xJ2=0.230; yJ2=0.150; zJ2=0.4865;
 xG2=0; yG2=0; zG2=0.230;
 m2=10; Ixx2=0.15; Iyy2=0.12; Izz2=0.05; Ixy2=0; Ixz2=0; Iyz2=0;   
 xJ3=0; yJ3=0; zJ3=0.475;
 xG3=-0.230; yG3=0; zG3=0;
 m3=5; Ixx3=0.05; Iyy3=0.05; Izz3=0.05; Ixy3=0; Ixz3=0; Iyz3=0;  
 xJ4=-0.230; yJ4=0.300; zJ4=0;
 xG4=0; yG4=0.450; zG4=0;
 m4=3; Ixx4=0.04; Iyy4=0.01; Izz4=0.04; Ixy4=0; Ixz4=0; Iyz4=0; 
 xJ5=0; yJ5=0.300; zJ5=0;
 xG5=0; yG5=0.300; zG5=0;
 m5=1; Ixx5=0.01; Iyy5=0.002; Izz5=0.01; Ixy5=0; Ixz5=0; Iyz5=0;
 xJ6=0; yJ6=0.065; zJ6=0;
 xG6=0; yG6=0.075; zG6=0;
 m6=1; Ixx6=0.005; Iyy6=0.005; Izz6=0.005; Ixy6=0; Ixz6=0; Iyz6=0;

% Inputs
g = 9.81; 
t1= 0; % Remember: it has been defined in the function !!!!!
t2= 0; 
t3= 0; 
t4= 0; 
t5= 0; 
t6= 0;

% ODE Solve
tic
z0 = rand(1,12); %[0 0 0 0 0 0 0 0 0 0 0 0].';
option = odeset('maxstep',0.1);
[t,z] = ode45(@IRB_1600,[0:0.01:1],z0,option);
toc

E_error = Energy_test(z);

q1 = z(:,1);   dq1 = z(:,7);                          
q2 = z(:,2);   dq2 = z(:,8);                          
q3 = z(:,3);   dq3 = z(:,9);                          
q4 = z(:,4);   dq4 = z(:,10);                         
q5 = z(:,5);   dq5 = z(:,11);                         
q6 = z(:,6);   dq6 = z(:,12);                         

figure
hold on
grid on
plot(t,q1,'linewidth',3)
plot(t,q2,'linewidth',3)
plot(t,q3,'linewidth',3)
plot(t,q4,'linewidth',3)
plot(t,q5,'linewidth',3)
plot(t,q6,'linewidth',3)
legend ('q1','q2','q3','q4','q5','q6')

figure
hold on
grid on
plot(t,dq1,'linewidth',3)
plot(t,dq2,'linewidth',3)
plot(t,dq3,'linewidth',3)
plot(t,dq4,'linewidth',3)
plot(t,dq5,'linewidth',3)
plot(t,dq6,'linewidth',3)
legend ('dq1','dq2','dq3','dq4','dq5','dq6')

figure
hold on
grid on
plot(t,E_error,'linewidth',3)
legend ('Energy Error')
