clear all
close all
clc

syms g t1 t2 t3 t4 t5 t6 real
syms q1 q2 q3 q4 q5 q6 real % six rotations
syms dq1 dq2 dq3 dq4 dq5 dq6 real % six angular velocities
q = [q1;q2;q3;q4;q5;q6];
dq = [dq1;dq2;dq3;dq4;dq5;dq6];

syms xJ1 yJ1 zJ1 real
syms xG1 yG1 zG1 m1 Ixx1 Iyy1 Izz1 Ixy1 Ixz1 Iyz1 real
syms xJ2 yJ2 zJ2 real
syms xG2 yG2 zG2 m2 Ixx2 Iyy2 Izz2 Ixy2 Ixz2 Iyz2 real
syms xJ3 yJ3 zJ3 real
syms xG3 yG3 zG3 m3 Ixx3 Iyy3 Izz3 Ixy3 Ixz3 Iyz3 real
syms xJ4 yJ4 zJ4 real
syms xG4 yG4 zG4 m4 Ixx4 Iyy4 Izz4 Ixy4 Ixz4 Iyz4 real
syms xJ5 yJ5 zJ5 real
syms xG5 yG5 zG5 m5 Ixx5 Iyy5 Izz5 Ixy5 Ixz5 Iyz5 real
syms xJ6 yJ6 zJ6 real
syms xG6 yG6 zG6 m6 Ixx6 Iyy6 Izz6 Ixy6 Ixz6 Iyz6 real

% Rotations
R1 = [cos(q1) sin(q1) 0 ; -sin(q1) cos(q1) 0 ; 0 0 1];
R2 = [1 0 0 ; 0 cos(q2) sin(q2) ; 0 -sin(q2) cos(q2)];
R3 = [1 0 0 ; 0 cos(q3) sin(q3) ; 0 -sin(q3) cos(q3)];
R4 = [cos(q4) 0 -sin(q4) ; 0 1 0 ; sin(q4) 0 cos(q4)];
R5 = [1 0 0 ; 0 cos(q5) sin(q5) ; 0 -sin(q5) cos(q5)];
R6 = [cos(q6) 0 -sin(q6) ; 0 1 0 ; sin(q6) 0 cos(q6)];

% Angular Moments
I1 = [Ixx1 -Ixy1 -Ixz1 ; -Ixy1 Iyy1 -Iyz1 ; -Ixz1 -Iyz1 Izz1];
I2 = [Ixx2 -Ixy2 -Ixz2 ; -Ixy2 Iyy2 -Iyz2 ; -Ixz2 -Iyz2 Izz2];
I3 = [Ixx3 -Ixy3 -Ixz3 ; -Ixy3 Iyy3 -Iyz3 ; -Ixz3 -Iyz3 Izz3];
I4 = [Ixx4 -Ixy4 -Ixz4 ; -Ixy4 Iyy4 -Iyz4 ; -Ixz4 -Iyz4 Izz4];
I5 = [Ixx5 -Ixy5 -Ixz5 ; -Ixy5 Iyy5 -Iyz5 ; -Ixz5 -Iyz5 Izz5];
I6 = [Ixx6 -Ixy6 -Ixz6 ; -Ixy6 Iyy6 -Iyz6 ; -Ixz6 -Iyz6 Izz6];

% Positions
r_J1 = [xJ1;yJ1;zJ1];
r_G1 = r_J1 + R1.'*[xG1;yG1;zG1];
r_J2 = r_J1 + R1.'*[xJ2;yJ2;zJ2];
r_G2 = r_J2 + R1.'*R2.'*[xG2;yG2;zG2];
r_J3 = r_J2 + R1.'*R2.'*[xJ3;yJ3;zJ3];
r_G3 = r_J3 + R1.'*R2.'*R3.'*[xG3;yG3;zG3];
r_J4 = r_J3 + R1.'*R2.'*R3.'*[xJ4;yJ4;zJ4];
r_G4 = r_J4 + R1.'*R2.'*R3.'*R4.'*[xG4;yG4;zG4];
r_J5 = r_J4 + R1.'*R2.'*R3.'*R4.'*[xJ5;yJ5;zJ5];
r_G5 = r_J5 + R1.'*R2.'*R3.'*R4.'*R5.'*[xG5;yG5;zG5];
r_J6 = r_J5 + R1.'*R2.'*R3.'*R4.'*R5.'*[xJ6;yJ6;zJ6];
r_G6 = r_J6 + R1.'*R2.'*R3.'*R4.'*R5.'*R6.'*[xG6;yG6;zG6];

% Velocities
V_G1 = jacobian(r_G1,q)*dq;
V_G2 = jacobian(r_G2,q)*dq;
V_G3 = jacobian(r_G3,q)*dq;
V_G4 = jacobian(r_G4,q)*dq;
V_G5 = jacobian(r_G5,q)*dq;
V_G6 = jacobian(r_G6,q)*dq;

% Angular Velocities
w1 = R1*[0;0;dq1];
w2 = R2*w1 + R2*[dq2;0;0];
w3 = R3*w2 + R3*[dq3;0;0];
w4 = R4*w3 + R4*[0;dq4;0];
w5 = R5*w4 + R5*[dq5;0;0];
w6 = R6*w5 + R6*[0;dq6;0];

% Generalized Forces
Power = [0 0 t1]*w1 + [t2 0 0]*(w2-w1) + [t3 0 0]*(w3-w2) + [0 t4 0]*(w4-w3) + [t5 0 0]*(w5-w4) + [0 t6 0]*(w6-w5);
Q = jacobian(Power,dq).';
% Kinetic Energy
T = ((0.5*m1*V_G1.'*V_G1)+(0.5*w1.'*I1*w1)) + ((0.5*m2*V_G2.'*V_G2)+(0.5*w2.'*I2*w2)) + ((0.5*m3*V_G3.'*V_G3)+(0.5*w3.'*I3*w3)) + ((0.5*m4*V_G4.'*V_G4)+(0.5*w4.'*I4*w4)) + ((0.5*m5*V_G5.'*V_G5)+(0.5*w5.'*I5*w5)) + ((0.5*m6*V_G6.'*V_G6)+(0.5*w6.'*I6*w6));

% Potential Energy
K = [0 0 1]*((m1*g*r_G1) + (m2*g*r_G2) + (m3*g*r_G3) + (m4*g*r_G4) + (m5*g*R1.'*r_G5) + (m6*g*r_G6));

% Lagrangian
L = T - K;
tic
Mass = jacobian(jacobian(L,dq).',dq);
Bias = (jacobian(jacobian(L,dq).',q)*dq - jacobian(L,q).');
Damp = 0.5*jacobian(Bias,dq);
KKKK = subs (Bias , dq , zeros(6,1));
toc

% Function Generation
fid=fopen('Required_Matrices.m','w');

fprintf(fid,'M11 = %s; \n', char((Mass(1,1))) );
fprintf(fid,'M12 = %s; \n', char((Mass(1,2))) );
fprintf(fid,'M13 = %s; \n', char((Mass(1,3))) );
fprintf(fid,'M14 = %s; \n', char((Mass(1,4))) );
fprintf(fid,'M15 = %s; \n', char((Mass(1,5))) );
fprintf(fid,'M16 = %s; \n', char((Mass(1,6))) );
fprintf(fid,'\n');

fprintf(fid,'M21 = %s; \n', char((Mass(2,1))) );
fprintf(fid,'M22 = %s; \n', char((Mass(2,2))) );
fprintf(fid,'M23 = %s; \n', char((Mass(2,3))) );
fprintf(fid,'M24 = %s; \n', char((Mass(2,4))) );
fprintf(fid,'M25 = %s; \n', char((Mass(2,5))) );
fprintf(fid,'M26 = %s; \n', char((Mass(2,6))) );
fprintf(fid,'\n');

fprintf(fid,'M31 = %s; \n', char((Mass(3,1))) );
fprintf(fid,'M32 = %s; \n', char((Mass(3,2))) );
fprintf(fid,'M33 = %s; \n', char((Mass(3,3))) );
fprintf(fid,'M34 = %s; \n', char((Mass(3,4))) );
fprintf(fid,'M35 = %s; \n', char((Mass(3,5))) );
fprintf(fid,'M36 = %s; \n', char((Mass(3,6))) );
fprintf(fid,'\n');

fprintf(fid,'M41 = %s; \n', char((Mass(4,1))) );
fprintf(fid,'M42 = %s; \n', char((Mass(4,2))) );
fprintf(fid,'M43 = %s; \n', char((Mass(4,3))) );
fprintf(fid,'M44 = %s; \n', char((Mass(4,4))) );
fprintf(fid,'M45 = %s; \n', char((Mass(4,5))) );
fprintf(fid,'M46 = %s; \n', char((Mass(4,6))) );
fprintf(fid,'\n');

fprintf(fid,'M51 = %s; \n', char((Mass(5,1))) );
fprintf(fid,'M52 = %s; \n', char((Mass(5,2))) );
fprintf(fid,'M53 = %s; \n', char((Mass(5,3))) );
fprintf(fid,'M54 = %s; \n', char((Mass(5,4))) );
fprintf(fid,'M55 = %s; \n', char((Mass(5,5))) );
fprintf(fid,'M56 = %s; \n', char((Mass(5,6))) );
fprintf(fid,'\n');

fprintf(fid,'M61 = %s; \n', char((Mass(6,1))) );
fprintf(fid,'M62 = %s; \n', char((Mass(6,2))) );
fprintf(fid,'M63 = %s; \n', char((Mass(6,3))) );
fprintf(fid,'M64 = %s; \n', char((Mass(6,4))) );
fprintf(fid,'M65 = %s; \n', char((Mass(6,5))) );
fprintf(fid,'M66 = %s; \n', char((Mass(6,6))) );
fprintf(fid,'\n');

fprintf(fid,'C11 = %s; \n', char((Damp(1,1))) );
fprintf(fid,'C12 = %s; \n', char((Damp(1,2))) );
fprintf(fid,'C13 = %s; \n', char((Damp(1,3))) );
fprintf(fid,'C14 = %s; \n', char((Damp(1,4))) );
fprintf(fid,'C15 = %s; \n', char((Damp(1,5))) );
fprintf(fid,'C16 = %s; \n', char((Damp(1,6))) );
fprintf(fid,'\n');

fprintf(fid,'C21 = %s; \n', char((Damp(2,1))) );
fprintf(fid,'C22 = %s; \n', char((Damp(2,2))) );
fprintf(fid,'C23 = %s; \n', char((Damp(2,3))) );
fprintf(fid,'C24 = %s; \n', char((Damp(2,4))) );
fprintf(fid,'C25 = %s; \n', char((Damp(2,5))) );
fprintf(fid,'C26 = %s; \n', char((Damp(2,6))) );
fprintf(fid,'\n');

fprintf(fid,'C31 = %s; \n', char((Damp(3,1))) );
fprintf(fid,'C32 = %s; \n', char((Damp(3,2))) );
fprintf(fid,'C33 = %s; \n', char((Damp(3,3))) );
fprintf(fid,'C34 = %s; \n', char((Damp(3,4))) );
fprintf(fid,'C35 = %s; \n', char((Damp(3,5))) );
fprintf(fid,'C36 = %s; \n', char((Damp(3,6))) );
fprintf(fid,'\n');

fprintf(fid,'C41 = %s; \n', char((Damp(4,1))) );
fprintf(fid,'C42 = %s; \n', char((Damp(4,2))) );
fprintf(fid,'C43 = %s; \n', char((Damp(4,3))) );
fprintf(fid,'C44 = %s; \n', char((Damp(4,4))) );
fprintf(fid,'C45 = %s; \n', char((Damp(4,5))) );
fprintf(fid,'C46 = %s; \n', char((Damp(4,6))) );
fprintf(fid,'\n');

fprintf(fid,'C51 = %s; \n', char((Damp(5,1))) );
fprintf(fid,'C52 = %s; \n', char((Damp(5,2))) );
fprintf(fid,'C53 = %s; \n', char((Damp(5,3))) );
fprintf(fid,'C54 = %s; \n', char((Damp(5,4))) );
fprintf(fid,'C55 = %s; \n', char((Damp(5,5))) );
fprintf(fid,'C56 = %s; \n', char((Damp(5,6))) );
fprintf(fid,'\n');

fprintf(fid,'C61 = %s; \n', char((Damp(6,1))) );
fprintf(fid,'C62 = %s; \n', char((Damp(6,2))) );
fprintf(fid,'C63 = %s; \n', char((Damp(6,3))) );
fprintf(fid,'C64 = %s; \n', char((Damp(6,4))) );
fprintf(fid,'C65 = %s; \n', char((Damp(6,5))) );
fprintf(fid,'C66 = %s; \n', char((Damp(6,6))) );
fprintf(fid,'\n');

fprintf(fid,'K1 = %s; \n', char((KKKK(1))) );
fprintf(fid,'K2 = %s; \n', char((KKKK(2))) );
fprintf(fid,'K3 = %s; \n', char((KKKK(3))) );
fprintf(fid,'K4 = %s; \n', char((KKKK(4))) );
fprintf(fid,'K5 = %s; \n', char((KKKK(5))) );
fprintf(fid,'K6 = %s; \n', char((KKKK(6))) );
fprintf(fid,'\n');

fprintf(fid,'M  = [M11 M12 M13 M14 M15 M16;                      \n');
fprintf(fid,'      M21 M22 M23 M24 M25 M26;                      \n');
fprintf(fid,'      M31 M32 M33 M34 M35 M36;                      \n');
fprintf(fid,'      M41 M42 M43 M44 M45 M46;                      \n');
fprintf(fid,'      M51 M52 M53 M54 M55 M56;                      \n');
fprintf(fid,'      M61 M62 M63 M64 M65 M66];                   \n\n');

fprintf(fid,'C  = [C11 C12 C13 C14 C15 C16;                      \n');
fprintf(fid,'      C21 C22 C23 C24 C25 C26;                      \n');
fprintf(fid,'      C31 C32 C33 C34 C35 C36;                      \n');
fprintf(fid,'      C41 C42 C43 C44 C45 C46;                      \n');
fprintf(fid,'      C51 C52 C53 C54 C55 C56;                      \n');
fprintf(fid,'      C61 C62 C63 C64 C65 C66];                   \n\n');

fprintf(fid,'K  = [K1;K2;K3;K4;K5;K6];                         \n\n');

fprintf(fid,'Q  = [t1;t2;t3;t4;t5;t6];                         \n\n');

fprintf(fid,'%% M*ddq + C*dq + K = Q                           \n\n');
  
fclose(fid);