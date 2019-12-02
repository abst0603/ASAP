%MAKE DATASET
%Head
clc
clear all
numpoint = 15000;
temp1 = 0.5*randn(numpoint,3)+[0.5,0,0];
temp2 = 0.8*randn(numpoint,3)+[2,0.5,0];
temp3 = 0.5*randn(numpoint,3)+[0,2,1];
temp4 = randn(numpoint,3)+[2,0,1];
temp5 = 0.8*randn(numpoint,3)+[1,1,0.5];
temp = [temp1;temp2;temp3;temp4;temp5];
temp = temp(temp(:,1)>0 & temp(:,1)<2 & temp(:,2)>0 & temp(:,2)<2 & temp(:,3)>0 & temp(:,3)<2,:);
plot3(temp(:,1),temp(:,2),temp(:,3),'r.');

clear temp1;clear temp2;clear temp3; clear temp4; clear temp5;
dataset = temp(pdist2(temp,[1,1,0])<1,:);
dataset = dataset(pdist2(dataset,[0.7,0.7,0.4])>0.3,:);
dataset1 = dataset(pdist2(dataset,[1.2,1.2,0.4])>0.35,:);
clear dataset; clear temp;
% plot3(dataset(:,1),dataset(:,2),dataset(:,3),'r.');
%two dimensional part
numpoint2 = 900;
temp1 = -1 * rand(5*numpoint2,1)+0.05;
temp4 = 0.1*randn(numpoint2,1)+0.3;
temp5 = 0.1*randn(numpoint2,1);
temp6 = 0.1*randn(numpoint2,1)+0.5;
temp7 = 0.1*randn(numpoint2,1)+0.8;
temp8 = 0.1*randn(numpoint2,1)+1.1;
temp9 = [temp4;temp5;temp6;temp7;temp8];
temp2 = 0.95*sin(temp9*pi*2)+1;
temp3 = 0.95*cos(temp9*pi*2)+1;
temp = [temp3,temp2,temp1];
dataset = temp(pdist2(temp,[0,1,-0.5])>0.2,:);
dataset = dataset(pdist2(dataset,[2,1,-0.5])>0.3,:);
dataset = dataset(pdist2(dataset,[1,2,-0.4])>0.2,:);
dataset = dataset(pdist2(dataset,[1,0,-1])>0.6,:);
dataset = [dataset;dataset1];
% plot3(dataset(:,1),dataset(:,2),dataset(:,3),'r.');
clear temp; clear temp1;clear temp2;clear temp3;clear temp4;clear temp5; clear temp6;clear temp7; clear temp8; clear temp9;
clear dataset1;
%tails
numpoint2 = 50;
temp1 = -1 * rand(numpoint2,1)-0.9;
temp4 = 0.008*rand(numpoint2,1)+0.3;
temp2 = 0.95*sin(temp4*pi*2)+1;
temp3 = 0.95*cos(temp4*pi*2)+1;
temp = [temp3,temp2,temp1];
dataset = [dataset;temp];
clear temp1; clear temp2; clear temp3; clear temp4; clear temp;
% plot3(dataset(:,1),dataset(:,2),dataset(:,3),'r.');
numpoint2 = 50;
temp1 = -1 * rand(numpoint2,1)-1;
temp4 = 0.01*rand(numpoint2,1)+0.75;
temp2 = 0.95*sin(temp4*pi*2)+1;
temp3 = 0.95*cos(temp4*pi*2)+1;
temp = [temp3,temp2,temp1];
dataset = [dataset;temp];
plot3(dataset(:,1),dataset(:,2),dataset(:,3),'r.');
clear temp1; clear temp2; clear temp3; clear temp4; clear temp;