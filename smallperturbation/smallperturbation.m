
%% Linearization around q1*=q2*=q3*=q4*=0 
%Xdot=A(x1....xn)

clc;clear;
syms k m l g M c S

A=[0 1 0 0 0 0 0 0
    -3*(g*l/2+k/m) -3*c/(m*l^2) 3*k/m 0 0 0 0 0
    0 0 0 1 0 0 0 0
    3*k/m 0 -3*(g*l/2+2*k/m) -3*c/(m*l^2) 3*k/m 0 0 0
    0 0 0 0 0 1 0 0
    0 0 3*k/m 0 -3*(g*l/2+2*k/m) -3*c/(m*l^2) 3*k/m 0
    0 0 0 0 0 0 0 1
    0 0 0 0 3*k/m 0 -3*(g*l/2+k/m) -3*c/(m*l^2)]

B=[0 0 0 0 0 0 0 M].'

C=eye(8)

D=0

eigenvluesA=eig(A)
polyA = charpoly(A,S)

% initial condition
data = {m, M, k, l, g, c};
datn = {0.5, 0, 0.75, 0.5, 9.81, 0.1};
A1 = subs(A, data, datn)
eig(A1)
B1 = subs(B, data, datn)
A1= [ 0     1         0     0         0     0         0     0
 -4743/400 -12/5       9/2     0         0     0         0     0
         0     0         0     1         0     0         0     0
       9/2     0 -6543/400 -12/5       9/2     0         0     0
         0     0         0     0         0     1         0     0
         0     0       9/2     0 -6543/400 -12/5       9/2     0
         0     0         0     0         0     0         0     1
         0     0         0     0       9/2     0 -4743/400 -12/5]
     
B1=[0 ;0 ;0 ;0; 0; 0; 0; 0]
     
eigenA1=eig(A1)
polyA1 = charpoly(A1,S)

%% Linearization around q1*=q2*=q3*=q4*=pi 
%Xdot=F(x1.....xn)

F=[0 1 0 0 0 0 0 0
    -3*(-g*l/2+k/m) -3*c/(m*l^2) 3*k/m 0 0 0 0 0
    0 0 0 1 0 0 0 0
    3*k/m 0 -3*(-g*l/2+2*k/m) -3*c/(m*l^2) 3*k/m 0 0 0
    0 0 0 0 0 1 0 0
    0 0 3*k/m 0 -3*(-g*l/2+2*k/m) -3*c/(m*l^2) 3*k/m 0
    0 0 0 0 0 0 0 1
    0 0 0 0 3*k/m 0 -3*(-g*l/2+k/m) -3*c/(m*l^2)]

eigenvluesF=eig(F)
polyF = charpoly(F,S)

% initial condition
data = {m, M, k, l, g, c};
datn = {0.5, 0, 0.75, 0.5, 9.81, 0.1};
F1 = subs(F, data, datn)
eig(F1)
F1= [ 0       1         0        0         0     0         0     0
 1143/400    -12/5       9/2     0         0     0         0     0
         0     0         0       1         0     0         0     0
       9/2     0 -657/400     -12/5       9/2      0         0     0
         0     0         0      0         0     1         0     0
         0     0       9/2      0   -657/400  -12/5       9/2     0
         0     0         0      0         0     0         0     1
         0     0         0      0       9/2     0 1143/400 -12/5]
    
     
eigenF1=eig(F1)
polyF1 = charpoly(F1,S)








