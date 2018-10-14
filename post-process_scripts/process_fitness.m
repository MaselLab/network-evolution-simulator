function fitness=process_fitness(sim_id,evo_steps,window_size)

%% load data
file_name=strcat('evo_summary_',num2str(sim_id),'.txt');
fid=fopen(file_name,'r');
data=textscan(fid,'%d %d %d %d %s %s %f %f %f %f %f %f %d %d %d %d\n',evo_steps,'HeaderLines',1);
fclose(fid);
close all

%% calcuate mean and variance
fitness=zeros(1,6);
fitness(1)=mean(data{8}(end-window_size+1:end)); 
fitness(2)=mean(data{9}(end-window_size+1:end));
fitness(3)=mean(data{7}(end-window_size+1:end));
fitness(4)=var(data{8}(end-window_size+1:end));
fitness(5)=var(data{9}(end-window_size+1:end));
fitness(6)=var(data{7}(end-window_size+1:end));

%% plot fitness vs evolutionary steps
plot(data{8},'b','DisplayName','env1 fitness'); 
hold on
plot(data{9},'r','DisplayName','env2 fitness');
plot(data{7},'m','DisplayName','mean fitness');
plot(evo_steps-window_size+1:evo_steps,data{7}(end-window_size+1:end),'k','DisplayName','mean f being evaluated');
legend('show','Location','best');
ax=gca;
ax.YLabel.String='fitness';
ax.XLabel.String='evoluationary steps';
saveas(gcf,'evolutionary_trajectory.jpg');
end