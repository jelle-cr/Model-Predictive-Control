close all;
clear all;
clc

Ts = 0.1;

A1 = 0;
B1 = 1;
C1 = 1;
D1 = 0;
sysc1 = ss(A1,B1,C1,D1);
sysd1 = c2d(sysc1,Ts);
A1 = sysd1.A;
B1 = sysd1.B;
C1 = sysd1.C;
D1 = sysd1.D;

Q1 = 1;
R1 = 1;
K1 = dlqr(A1,B1,Q1,R1);

A2 = [0 1; 0 0];
B2 = [0; -1];
C2 = [1 0; 0 1];
D2 = 0;
sysc2 = ss(A2,B2,C2,D2);
sysd2 = c2d(sysc2,Ts);
A2 = sysd2.A;
B2 = sysd2.B;
C2 = sysd2.C;
D2 = sysd2.D;
%y2 = lsim(sysd2,

Q2 = eye(2);
R2 = 1;
K2 = dlqr(A2,B2,Q2,R2);