clear all;
close all;
clc;

N = 4;

%% First order model
load modelFirstOrder.mat
[n,m] = size(B);

%% Question 2.2a
umin = -0.3;
umax = 0.3;
xmin = -10^9;
xmax = 10^9;

%[CcalTest, DcalTest, EcalTest, McalTest] = caligraphicMatricesTest(umin,umax,xmin,xmax,N,n,m);
[Ccal, Dcal, Ecal, Mcal] = caligraphicMatrices(umin,umax,xmin,xmax,N,n,m);

