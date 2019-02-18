function [means,vars]=summarize_parameters(ids,N_samples)
%%
%This function calcualted the mean and variance of the sampled parameters
%of network motifs.
%%
%number of simulation replicates
n=length(ids);

%temporary storage for sampled parameters. Each row is samples taken from
%one simulation replicate
sinal_kd=zeros(n,N_samples); %kd of TF
effector_l=zeros(n,N_samples); %locus length
effector_md=zeros(n,N_samples); %mRNA degradation rate
effector_pd=zeros(n,N_samples); %protein degradtion rate
effector_tl=zeros(n,N_samples); %protein synthesis rate
effector_a2i=zeros(n,N_samples); %transition rate of active promoter state to intermediate promoter state 
slow_TF_l=zeros(n,N_samples);
slow_TF_md=zeros(n,N_samples);
slow_TF_pd=zeros(n,N_samples);
slow_TF_tl=zeros(n,N_samples);
slow_TF_a2i=zeros(n,N_samples);
slow_TF_kd=zeros(n,N_samples);
fast_TF_l=zeros(n,N_samples);
fast_TF_md=zeros(n,N_samples);
fast_TF_pd=zeros(n,N_samples);
fast_TF_tl=zeros(n,N_samples);
fast_TF_a2i=zeros(n,N_samples);
fast_TF_kd=zeros(n,N_samples);

%loop through simulation replicates
for i=1:n    
    %cd into a replicate
    cd(strcat('rep',num2str(ids(i))));
    cd('result');
    %load result
    data=load('parameters.txt','r'); 
    
    %sort samples of each parameter
    for j=1:N_samples
        sinal_kd(i,j)=data(4*(j-1)+1,8);
        effector_a2i(i,j)=data(4*(j-1)+2,3);
        effector_md(i,j)=data(4*(j-1)+2,4);
        effector_tl(i,j)=data(4*(j-1)+2,5);
        effector_pd(i,j)=data(4*(j-1)+2,6);
        effector_l(i,j)=data(4*(j-1)+2,7);
        fast_TF_a2i(i,j)=data(4*(j-1)+3,3);        
        fast_TF_md(i,j)=data(4*(j-1)+3,4); 
        fast_TF_tl(i,j)=data(4*(j-1)+3,5); 
        fast_TF_pd(i,j)=data(4*(j-1)+3,6); 
        fast_TF_l(i,j)=data(4*(j-1)+3,7); 
        fast_TF_kd(i,j)=data(4*(j-1)+3,8); 
        slow_TF_a2i(i,j)=data(4*j,3);        
        slow_TF_md(i,j)=data(4*j,4); 
        slow_TF_tl(i,j)=data(4*j,5); 
        slow_TF_pd(i,j)=data(4*j,6); 
        slow_TF_l(i,j)=data(4*j,7); 
        slow_TF_kd(i,j)=data(4*j,8);
    end 
    cd ..;
    cd ..;
end

%average parameters over one batch of sample taken from
%each simulation replicates (average of a column), then average of all batches of samples.
mean_signal=mean(mean(sinal_kd,1));
mean_effector(1)=mean(mean(effector_a2i,1));
mean_effector(2)=mean(mean(effector_md,1));
mean_effector(3)=mean(mean(effector_tl,1));
mean_effector(4)=mean(mean(effector_pd,1));
mean_effector(5)=mean(mean(effector_l,1));
mean_slow_TF(1)=mean(mean(slow_TF_a2i,1));
mean_slow_TF(2)=mean(mean(slow_TF_md,1));
mean_slow_TF(3)=mean(mean(slow_TF_tl,1));
mean_slow_TF(4)=mean(mean(slow_TF_pd,1));
mean_slow_TF(5)=mean(mean(slow_TF_l,1));
mean_slow_TF(6)=mean(mean(slow_TF_kd,1));
mean_fast_TF(1)=mean(mean(fast_TF_a2i,1));
mean_fast_TF(2)=mean(mean(fast_TF_md,1));
mean_fast_TF(3)=mean(mean(fast_TF_tl,1));
mean_fast_TF(4)=mean(mean(fast_TF_pd,1));
mean_fast_TF(5)=mean(mean(fast_TF_l,1));
mean_fast_TF(6)=mean(mean(fast_TF_kd,1));

means=[mean_fast_TF;mean_slow_TF;mean_effector, NAN; NAN,NAN,NAN,NAN,NAN,mean_signal];

%calculate the variance of parameters over one batch of sample taken from
%each simulation replicates (average of a column), then average of all batches of samples.
var_signal=mean(var(sinal_kd,0,1));
var_effector(1)=mean(var(effector_a2i,0,1));
var_effector(2)=mean(var(effector_md,0,1));
var_effector(3)=mean(var(effector_tl,0,1));
var_effector(4)=mean(var(effector_pd,0,1));
var_effector(5)=mean(var(effector_l,0,1)./(mean(effector_l,1)).^2); %for locus length we calculate the coefficient of variance instead of variance
var_slow_TF(1)=mean(var(slow_TF_a2i,0,1));
var_slow_TF(2)=mean(var(slow_TF_md,0,1));
var_slow_TF(3)=mean(var(slow_TF_tl,0,1));
var_slow_TF(4)=mean(var(slow_TF_pd,0,1));
var_slow_TF(5)=mean(var(slow_TF_l,0,1)./(mean(slow_TF_l,1)).^2);
var_slow_TF(6)=mean(var(slow_TF_kd,0,1));
var_fast_TF(1)=mean(var(fast_TF_a2i,0,1));
var_fast_TF(2)=mean(var(fast_TF_md,0,1));
var_fast_TF(3)=mean(var(fast_TF_tl,0,1));
var_fast_TF(4)=mean(var(fast_TF_pd,0,1));
var_fast_TF(5)=mean(var(fast_TF_l,0,1)./(mean(fast_TF_l,1)).^2);
var_fast_TF(6)=mean(var(fast_TF_kd,0,1));

vars=[var_fast_TF;var_slow_TF;var_effector, NAN; NAN,NAN,NAN,NAN,NAN,var_signal];