function [occurrence,contribution_to_near_AND]=calc_near_AND_gated_motifs(window_size)
%%
%The function calculate the proportions of evolutionary steps that contain
%AND-gated motifs and/or near-AND-gated motifs. And also indicate the
%proportions of near-AND-gated motifs that are fast-TF-controlled,
%slow-TF-controlled, or OR-gated.
%%

occurrence=zeros(1,16);
contribution_to_near_AND=zeros(1,6);


%load the number of near-AND-gated motifs in each evolutionary steps.    
temp1=load('N_near_AND_gate_motifs.txt');
%only count the occurrence in the last 10000 steps
temp1(1:end-window_size,:)=[];  

%calculate the total number of near-AND-gated isolated C1-ffls,
%FFL-in-diamonds, and isolated diamonds
near_and_gated=[sum(temp1(:,2:4),2),sum(temp1(:,6:8),2),sum(temp1(:,10:12),2)];    
near_AND_c1ffl=[sum(temp1(:,2:4),2),temp1(:,2:4)];
near_AND_ffl_in_diamond=[sum(temp1(:,6:8),2),temp1(:,6:8)];
near_AND_diamond=[sum(temp1(:,10:12),2),temp1(:,10:12)];

%find the steps that contain at least on near-AND-gated motif
near_AND_c1ffl(near_AND_c1ffl(:,1)==0,:)=[];
near_AND_ffl_in_diamond(near_AND_ffl_in_diamond(:,1)==0,:)=[];
near_AND_diamond(near_AND_diamond(:,1)==0,:)=[];

%calculate the proportions of near-AND-gated motifs that are fast-TF-controlled, slow-TF-controlled, or OR-gated.
contribution_to_near_AND(1:3)=mean([near_AND_c1ffl(:,2)./near_AND_c1ffl(:,1),near_AND_c1ffl(:,3)./near_AND_c1ffl(:,1),near_AND_c1ffl(:,4)./near_AND_c1ffl(:,1)],1);
contribution_to_near_AND(4:6)=mean([near_AND_ffl_in_diamond(:,2)./near_AND_ffl_in_diamond(:,1),near_AND_ffl_in_diamond(:,3)./near_AND_ffl_in_diamond(:,1),near_AND_ffl_in_diamond(:,4)./near_AND_ffl_in_diamond(:,1)],1);
contribution_to_near_AND(7:9)=mean([near_AND_diamond(:,2)./near_AND_diamond(:,1),near_AND_diamond(:,3)./near_AND_diamond(:,1),near_AND_diamond(:,4)./near_AND_diamond(:,1)],1);

%col1, 5, 9 of N_near_AND_gate_motifs.txt are the number of AND-gated
%motifs that contain only strong TFBSs in the effector gene. Call these
%AND-gated motifs strong-AND-gated 
strong_and=temp1(:,[1,5,9]);

%load the number of AND-gated motifs in each evolutionary step
temp2=load('N_motifs.txt');
temp2(1:end-window_size,:)=[];
and_gated=temp2(:,[15,24,33]);

%Call AND-gated motifs that depend on weak TFBSs in the effector gene
%weak-AND-gated motifs      
weak_and=and_gated-strong_and;

%count the steps that contain AND-gated alone, near AND-gated alone, and co-exist  
near_and_gated=near_and_gated>0;
and_gated=and_gated>0;
weak_and=weak_and>0;    

coexist_and_near_and=(near_and_gated+and_gated)==2;
coexist_and_near_and=sum(coexist_and_near_and,1)/window_size; 

either_and_near_and=(near_and_gated+and_gated)~=0;    
%also count the steps that contain at least one of AND-gated isolated
%C1-FFLs, near-AND-gated isolated C1-FFLs, AND-gated isolated diamonds,
%and near-AND-gated isolated diamonds
either_c1ffl_diamond=(either_and_near_and(:,1)+either_and_near_and(:,3))~=0; 
either_and_near_and=sum(either_and_near_and,1)/window_size;
either_c1ffl_diamond=sum(either_c1ffl_diamond,1)/window_size;

natural_and_alone=(and_gated-near_and_gated)==1;
natural_and_alone=sum(natural_and_alone,1)/window_size;

near_and_alone=(near_and_gated-and_gated)==1;
near_and_alone=sum(near_and_alone,1)/window_size;

weak_and=sum(weak_and,1)/window_size;    

%pool results
occurrence(1:5)=[natural_and_alone(1),near_and_alone(1),coexist_and_near_and(1),either_and_near_and(1),weak_and(1)]; %about isoalted C1-FFLs
occurrence(6:10)=[natural_and_alone(2),near_and_alone(2),coexist_and_near_and(2),either_and_near_and(2),weak_and(2)]; %about FFL-in-diamonds
occurrence(11:15)=[natural_and_alone(3),near_and_alone(3),coexist_and_near_and(3),either_and_near_and(3),weak_and(3)]; %about isolated diamonds
occurrence(16)=either_c1ffl_diamond;
    
