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

%% Lagrangian
Lagr=T-U
fprintf('Lagr = \n'); pretty(simplify(Lagr)); 
fprintf('\n')

%% generalized momenta
syms pq1 pq2 pq3 pq4
Pq1=deriv(T,diff(q1,t))
Pq2=deriv(T,diff(q2,t))
Pq3=deriv(T,diff(q3,t))
Pq4=deriv(T,diff(q4,t))

A=[m*L^2/3 0 0 0;0 m*L^2/3 0 0;0 0 m*L^2/3 0;0 0 0 m*L^2/3]  % pi=A*dqi ==> dqi=inv(A)pi

dq1=[1 0 0 0]*inv(A)*[pq1;pq2;pq3;pq4]
dq2=[0 1 0 0]*inv(A)*[pq1;pq2;pq3;pq4]
dq3=[0 0 1 0]*inv(A)*[pq1;pq2;pq3;pq4]
dq4=[0 0 0 1]*inv(A)*[pq1;pq2;pq3;pq4]

%% hamiltonian function
H=pq1*wA+pq2*wB+pq3*wC+pq4*wD-Lagr;

H=subs(H,diff(q1,t),dq1);
H=subs(H,diff(q2,t),dq2);
H=subs(H,diff(q3,t),dq3);
H=subs(H,diff(q4,t),dq4);

fprintf('H = \n'); pretty(simplify(H)); 
fprintf('\n')

%% generalized noncoservative force Qn1
Qn1 = -3*c/(m*L^2)*pq1;
%% generalized noncoservative force Qn2
Qn2 = -3*c/(m*L^2)*pq2;
%% generalized noncoservative force Qn3
Qn3 = -3*c/(m*L^2)*pq3;
%% generalized noncoservative force Qn4
Qn4 = -3*c/(m*L^2)*pq4+M;

fprintf('Qn1 = \n'); pretty(simplify(Qn1)); 
fprintf('\n')
fprintf('Qn2 = \n'); pretty(simplify(Qn2)); 
fprintf('\n')
fprintf('Qn3 = \n'); pretty(simplify(Qn3)); 
fprintf('\n')
fprintf('Qn4 = \n'); pretty(simplify(Qn4)); 
fprintf('\n')


%% Canonical Hamilton Equation Of Motion 

fprintf( ' \n  \n Canonical Hamilton Equation Of Motion  \n'); 
fprintf('\n')

%first equation
dq1=deriv(H,pq1);
dq1 = simplify(dq1)

%second equation
dq1=deriv(H,pq2);
dq2 = simplify(dq2)

%third equation
dq3=deriv(H,pq3);
dq3 = simplify(dq3)


%fourth equation
dq4=deriv(H,pq4);
dq4 = simplify(dq4)

%fifth equation
dp1=-deriv(H,q1)+Qn1;
dp1 = simplify(dp1)

%sixth equation
dp2=-deriv(H,q2)+Qn2;
dp2 = simplify(dp2)

%seventh equation
dp3=-deriv(H,q3)+Qn3;
dp3 = simplify(dp3)

%eighth equation
dp4=-deriv(H,q4)+Qn4;
dp4 = simplify(dp4)




