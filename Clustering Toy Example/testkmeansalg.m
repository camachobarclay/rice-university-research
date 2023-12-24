clear all;
close all;
clc;



cd(matlabroot)
cd ..
cd ..
cd 'E:\Documents\Rice University\Research\Intro_to_Cluster_Analysis'
%cd ~/documents/Rice' University'/Research/Intro' to Cluster Analysis'/


X = importdata('s1.txt')/10000;
tol =10^(-12);
maxiter = 1000;
mu1 = [24, 84;
14,56;
17,35;
42,79;
34,56;
40,40;
32,17;
67,86;
61,56;
63,40;
51,18;
83,73;
86,55;
80,33;
85,15];

% X = [randn(100,2)+ones(100,2);...
%      randn(100,2)-ones(100,2)];
% opts = statset('Display','final');
k = 15;

n = 4;
timevalues = zeros(2*n,1);

 
colorval = zeros(1,3);

hold on
for t=1:n
    tic;
[IDX, C] = kmeans(X,k,'distance','cityblock');
timeval = toc;
timevalues(t) = timeval;
[C2, mu2, MU, superrval, timeval] = my_kmeans(X,k,maxiter, tol);
timevalues(t+4) = timeval;
    for j=1:k
        figure(2*t-1)
        hold on
        colorval = rand(1,3);
        plot(X(IDX==j,1),X(IDX==j,2),'o', 'color',...
            colorval,'MarkerSize', 5)
        plot(X(IDX==j,1),X(IDX==j,2),'.', 'color',...
            colorval,'MarkerSize', 2)
        plot(C(j,1),C(j,2),'kx', 'MarkerSize',12,'LineWidth',2)
        plot(C(j,1),C(j,2),'o', 'color',...
            colorval, 'MarkerSize',12,'LineWidth',2)
        
        
        figure(2*t)
        hold on
        colorval = rand(1,3);
        colorval2 = rand(1,3);
        plot(X(C2==j,1),X(C2==j,2),'o', 'color',...
            colorval,'MarkerSize', 5)
        plot(X(C2==j,1),X(C2==j,2),'.', 'color',...
            colorval, 'MarkerSize', 2)
        plot(mu2(j,1),mu2(j,2),'kx', 'MarkerSize',12,'LineWidth',2)
        plot(mu2(j,1),mu2(j,2),'o', 'color',...
            colorval, 'MarkerSize',12,'LineWidth',2)

        
    end
        figure(2*t-1)
        title(['MATLAB''s kmeans ell-1 Clustering of the Data,'...
            'Test Run #' num2str(t)])
        xlabel('x-axis')
        ylabel('y-axis')
        
        figure(2*t)
        title(['User Defined kmeans ell-1 Clustering of the Data,'...
            'Test Run #' num2str(t)])
        xlabel('x-axis')
        ylabel('y-axis')
        figure(2*n+1)
        hold on
        colorval = rand(1,3);
        plot(0:length(superrval)-1,superrval, 'color',...
            colorval, 'linewidth', 3);
end

figure(2*n+1)
title(['Maxiumum Change in Magnitude Between Iterations',...
    'of the User''s kmeans ell-1 Algorithm'])
xlabel('Number of Iterations')
ylabel('Value')
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4')
