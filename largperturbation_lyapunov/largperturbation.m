%% Advanced dyanmic project
%M=0
%teta4=arbitrary

clear; 
clc;
close all,

syms L m k c g t M  

q1 = sym('q1(t)');
q2 = sym('q2(t)');
q3 = sym('q3(t)');
q4 = sym('q4(t)');

rA= [L*sin(q1)/2,L*cos(q1)/2,0];
rB= [L*sin(q2)/2,L*cos(q2)/2,0];
rC= [L*sin(q3)/2,L*cos(q3)/2,0];
rD= [L*sin(q4)/2,L*cos(q4)/2,0];

wA=diff(q1,t)
wB=diff(q2,t)
wC=diff(q3,t)
wD=diff(q4,t)

vA=diff(rA,t)
vB=diff(rB,t)
vC=diff(rC,t)
vD=diff(rD,t)

%% kinetic energy of 1   
%where I is around center of rotation Io=m*L^2/3

T1 = m*L^2*wA*wA.'/6;

%% kinetic energy of 2
T2 = m*L^2*wB*wB.'/6;

%% kinetic energy of 3
T3 = m*L^2*wC*wC.'/6;

%% kinetic energy of 4
T4 = m*L^2*wD*wD.'/6;

fprintf('T1 = %s \n',char(T1))
fprintf('T2 = %s \n',char(T2))
fprintf('T3 = %s \n',char(T3))
fprintf('T4 = %s \n',char(T4))



%% total kinetic energy
T = expand(T1 + T2 + T3 + T4);
fprintf('T = \n'); pretty(simplify(T)); 
fprintf('\n')

%% total potential energy
U1=m*g*L*(1-cos(q1))/2
U2=m*g*L*(1-cos(q2))/2
U3=m*g*L*(1-cos(q3))/2
U4=m*g*L*(1-cos(q4))/2
Ue=k*L^2*((q1-q2)^2+(q2-q3)^2+(q3-q4)^2)/2
U=U1+U2+U3+U4+Ue
fprintf('U = \n'); pretty(simplify(U)); 
fprintf('\n')

%dU/dq
Uq1 = deriv(U, q1);
Uq2 = deriv(U, q2);
Uq3 = deriv(U, q3);
Uq4 = deriv(U, q4);
fprintf('dU/dq1 = %s \n',char(Uq1))
fprintf('\n')
fprintf('dU/dq2 = %s \n',char(Uq2)) 
fprintf('\n')
fprintf('dU/dq3 = %s \n',char(Uq3))
fprintf('\n')
fprintf('dU/dq4 = %s \n',char(Uq4)) 
fprintf('\n')



% dT/d(dq)
Tdq1 = deriv(T, diff(q1,t));
Tdq2 = deriv(T, diff(q2,t));
Tdq3 = deriv(T, diff(q3,t));
Tdq4 = deriv(T, diff(q4,t));
fprintf('dT/d(dq1) = \n'); pretty(simplify(Tdq1));
fprintf('\n')
fprintf('dT/d(dq2) = \n'); pretty(simplify(Tdq2)); 
fprintf('\n')
fprintf('dT/d(dq3) = \n'); pretty(simplify(Tdq3));
fprintf('\n')
fprintf('dT/d(dq4) = \n'); pretty(simplify(Tdq4)); 
fprintf('\n')

% d(dT/d(dq))/dt
Tt1 = diff(Tdq1, t);
Tt2 = diff(Tdq2, t);
Tt3 = diff(Tdq3, t);
Tt4 = diff(Tdq4, t);
fprintf('d dT/d(dq1)/dt = \n'); 
pretty(simplify(Tt1)); 
fprintf('\n')
fprintf('d dT/d(dq2)/dt = \n'); 
pretty(simplify(Tt2)); 
fprintf('\n')
fprintf('d dT/d(dq3)/dt = \n'); 
pretty(simplify(Tt3)); 
fprintf('\n')
fprintf('d dT/d(dq4)/dt = \n'); 
pretty(simplify(Tt4)); 
fprintf('\n')

% dT/dq
Tq1 = deriv(T, q1);
Tq2 = deriv(T, q2);
Tq3 = deriv(T, q3);
Tq4 = deriv(T, q4);
fprintf('dT/dq1 = %s \n',char(Tq1))
fprintf('\n')
fprintf('dT/dq2 = %s \n',char(Tq2)) 
fprintf('\n')
fprintf('dT/dq3 = %s \n',char(Tq3))
fprintf('\n')
fprintf('dT/dq4 = %s \n',char(Tq4)) 
fprintf('\n')

%Raileigh's dissipation function 
F1= c*wA*wA.'/2;
F2= c*wB*wB.'/2;
F3= c*wC*wC.'/2;
F4= c*wD*wD.'/2;

%viscus friction generalized forces
Fdq1 = deriv(F1, diff(q1,t));
Fdq2 = deriv(F2, diff(q2,t));
Fdq3 = deriv(F3, diff(q3,t));
Fdq4 = deriv(F4, diff(q4,t));

%% left hand side of Lagrange's eom
LHS1 = Tt1 - Tq1 + Fdq1;
LHS2 = Tt2 - Tq2 + Fdq2;
LHS3 = Tt3 - Tq3 + Fdq3;
LHS4 = Tt4 - Tq4 + Fdq4;

%% generalized noncoservative active forces

G1 = [0 0 0]; 
G2 = [0 0 0];
G3 = [0 0 0]; 
G4 = [2*M/(L*cos(q4)) 0  0];
%% partial derivatives
rA_1 = deriv(rA, q1); rB_1 = deriv(rB, q1);rC_1 = deriv(rC, q1); rD_1 = deriv(rD, q1);
rA_2 = deriv(rA, q2); rB_2 = deriv(rB, q2);rC_2 = deriv(rC, q2); rD_2 = deriv(rD, q2);
rA_3 = deriv(rA, q3); rB_3 = deriv(rB, q3);rC_3 = deriv(rC, q3); rD_3 = deriv(rD, q3);
rA_4 = deriv(rA, q4); rB_4 = deriv(rB, q4);rC_4 = deriv(rC, q4); rD_4 = deriv(rD, q4);
%% generalized active force Q1
Q1 = rA_1*G1.'+rB_1*G2.'+rC_1*G3.'+rD_1*G4.'-Uq1;
%% generalized active force Q2
Q2 = rA_2*G1.'+rB_2*G2.'+rC_2*G3.'+rD_2*G4.'-Uq2;
%% generalized active force Q3
Q3 = rA_3*G1.'+rB_3*G2.'+rC_3*G3.'+rD_3*G4.'-Uq3;
%% generalized active force Q4
Q4 = rA_4*G1.'+rB_4*G2.'+rC_4*G3.'+rD_4*G4.'-Uq4;
fprintf('Q1 = \n'); pretty(simplify(Q1)); 
fprintf('\n')
fprintf('Q2 = \n'); pretty(simplify(Q2)); 
fprintf('\n')
fprintf('Q3 = \n'); pretty(simplify(Q3)); 
fprintf('\n')
fprintf('Q4 = \n'); pretty(simplify(Q4)); 
fprintf('\n')
%% first Lagrange's equation of motion
Lagrange1 = LHS1-Q1
%% second Lagrange's equation of motion
Lagrange2 = LHS2-Q2
%% third Lagrange's equation of motion
Lagrange3 = LHS3-Q3
%% fourth Lagrange's equation of motion
Lagrange4 = LHS4-Q4

ql = {diff(q1,t,2), diff(q2,t,2),diff(q3,t,2),diff(q4,t,2)...
    diff(q1,t), diff(q2,t),diff(q3,t),diff(q4,t),q1, q2,q3,q4};    
qf = ...
{'ddq1', 'ddq2','ddq3','ddq4' 'y(2)', 'y(4)','y(6)','y(8)', 'y(1)', 'y(3)', 'y(5)','y(7)'};

Lagra1 = subs(Lagrange1, ql, qf);
Lagra2 = subs(Lagrange2, ql, qf);
Lagra3 = subs(Lagrange3, ql, qf);
Lagra4 = subs(Lagrange4, ql, qf);

%% solve e.o.m. for ddq1, ddq2 , ddq3 , ddq4
sol = solve(Lagra1,Lagra2,Lagra3,Lagra4,'ddq1, ddq2, ddq3 , ddq4');
dy2 = sol.ddq1
dy4 = sol.ddq2
dy6 = sol.ddq3
dy8 = sol.ddq4

%% lyapunov Function
lyapunov1=T+U;
fprintf('lyapunov1 = \n'); pretty(simplify(lyapunov1)); 
fprintf('\n')

vdot=deriv(lyapunov1,q1)*diff(q1,t)+deriv(lyapunov1,q2)*diff(q2,t)+deriv(lyapunov1,q3)*diff(q3,t)+deriv(lyapunov1,q4)*diff(q4,t)+deriv(lyapunov1,diff(q1,t))*diff(q1,t,2)+deriv(lyapunov1,diff(q2,t))*diff(q2,t,2)+deriv(lyapunov1,diff(q3,t))*diff(q3,t,2)+deriv(lyapunov1,diff(q4,t))*diff(q4,t,2);

q2 = {diff(q1,t,2), diff(q2,t,2),diff(q3,t,2),diff(q4,t,2)...
    diff(q1,t), diff(q2,t),diff(q3,t),diff(q4,t),q1, q2,q3,q4};    
qf2 = ...
{dy2,dy4 ,dy6,dy8, 'y(2)', 'y(4)','y(6)','y(8)', 'y(1)', 'y(3)', 'y(5)','y(7)'};

V = subs(lyapunov1, q2, qf2);
fprintf('V = \n'); pretty(simplify(V)); 
fprintf('\n')
vdot=subs(vdot, q2, qf2);
data = {M};
datn = {0};
vdot = subs(vdot, data, datn);
fprintf('vdot = \n'); pretty(simplify(vdot)); 
fprintf('\n')


%matrix of Vdot with initialdata = {m=0.5, M=0, k=0.75, L=0.5, g=9.81, c=0.1};

Q=[0 0 0 0 0 0 0 0
0 .1 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 .1 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 .1 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 .1]

A1= [ 0     1         0     0         0     0         0     0
 -4743/400 -12/5       9/2     0         0     0         0     0
         0     0         0     1         0     0         0     0
       9/2     0 -6543/400 -12/5       9/2     0         0     0
         0     0         0     0         0     1         0     0
         0     0       9/2     0 -6543/400 -12/5       9/2     0
         0     0         0     0         0     0         0     1
         0     0         0     0       9/2     0 -4743/400 -12/5]
     
     
     
P=lyap(A1,Q);
P=vpa(P)
eigenvaluesP=eig(P)
%P is positive difinite hence the system is asym.stable



%% calculate V for inital condition
syms y1 y2 y3 y4 y5 y6 y7 y8

VA1=[y1 y2 y3 y4 y5 y6 y7 y8]*P*[y1 y2 y3 y4 y5 y6 y7 y8].';
VA1=simplify(VA1)

VdotA1=-[y1 y2 y3 y4 y5 y6 y7 y8]*Q*[y1 y2 y3 y4 y5 y6 y7 y8].';
VdotA1=simplify(VdotA1)


%% plot V in domain of attraction
m=0.5;
L=0.5;
g=9.81;


%state 1
data = {y2, y4, y6, y7, y8 };
datn = {0, 0, 0 , 0, 0};
Vs1 = subs(VA1, data, datn);
region1=Vs1-4*m*g*L == 0
figure
f1 =region1;
interval = [-100 100 -100 100 -100 100];
fimplicit3(f1,interval)

%state 2
data = {y2, y4, y5, y6, y8 };
datn = {0, 0, 0 , 0, 0};
Vs2 = subs(VA1, data, datn);
region2=Vs2-4*m*g*L == 0
figure
f2 =region2;
interval = [-100 100 -100 100 -100 100];
fimplicit3(f2,interval)

%state3
data = {y2, y3 ,y4 ,y6, y8 };
datn = {0, 0, 0 , 0, 0};
Vs3 = subs(VA1, data, datn);
region3=Vs3-4*m*g*L == 0
figure
f3 =region3;
interval = [-100 100 -100 100 -100 100];
fimplicit3(f3,interval)

%state4
data = {y1 ,y2, y4 ,y6, y8 };
datn = {0, 0, 0 , 0, 0};
Vs4 = subs(VA1, data, datn);
region4=Vs4-4*m*g*L == 0
figure
f4 =region4;
interval = [-100 100 -100 100 -100 100];
fimplicit3(f4,interval)





