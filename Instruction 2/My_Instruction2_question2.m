clear all;
close all;
clc;

N = 4;

%% First order model
load modelSecondOrder.mat
[n,m] = size(B);

%% Question 2.2a
%[CcalTest, DcalTest, EcalTest, McalTest] = caligraphicMatricesTest(umin,umax,xmin,xmax,N,n,m);
[Ccal, Dcal, Ecal, Mcal] = caligraphicMatrices(umin,umax,xmin,xmax,N,n,m);

%% Question 2.2b
[phi, gamma] = predictionModel(A,B,N,n,m);
L = Mcal*gamma + Ecal;
W = -Dcal-Mcal*phi;

