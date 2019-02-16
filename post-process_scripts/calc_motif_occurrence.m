function motif_occurrence=calc_motif_occurrence(window_size)  
%%load data
data=load('N_motifs.txt'); 
motif_numbers=data(end-window_size+1:end,1:36); 
%%check for presence of a motif
motif_presence=motif_numbers>0;    
%%calculate occurrence
motif_occurrence=sum(motif_presence,1)./window_size;  
end