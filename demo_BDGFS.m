addpath(genpath(pwd));

clear;
clc;
load('./datasets/lungd.mat');


NewFeaNum=[20 30 40 50 60 70 80 90 100];
nc=length(unique(gnd));
alpha = 1e3;%block-diagonal
gamma = 1e-2;%equal of U-Z
lambda = 1;%sparsity
beta = 0.001;%0.001;%graph penalty

pars=[1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3];
% pars=[2e-1,4e-1,6e-1,8e-1,1e0,2e0,4e0,6e0,8e0];
tic;


for i1=1:length(pars)
    alpha=pars(i1);
    for i2=1:length(pars)
        gamma=pars(i2);
        
        fprintf('alpha=%.4f, gamma=%.4f\n',alpha,gamma);
        [W,obj]=BDGFS(fea,alpha,gamma,lambda,beta,nc);
        toc;
        [~, idx] = sort(sum(W.*W,2),'descend');
        for i=1:size(NewFeaNum,2) 
            Newfea=fea(:,idx(1:NewFeaNum(i)));
            label=litekmeans(Newfea,nc,'Replicates',20);
            [Acc(i),Nmi(i),~]=ClusteringMeasure(gnd,label);
        end
    end
end

