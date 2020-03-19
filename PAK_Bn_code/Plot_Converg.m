% Plot convergence history

iter_mat = cell2mat(iter_cell);
iter_sum = [0; cumsum(iter_mat)];
chg_pt = cumsum(iter_mat)+N_iniDOE;
total_iter = sum(iter_mat);
figure()
for i=1:k 
    plot(N_iniDOE+iter_sum(i):N_iniDOE+iter_sum(i+1),P_f_cell{i,1}/P_f_MCS,'LineWidth',1.5,'Color','m')
    if i==1
        hold on
    end   
    if i==1
        plot(N_iniDOE:N_iniDOE+total_iter,ones(1,total_iter+1),'LineWidth',1.5,'Color','k') % "Exact" sol (MCS)
    end 
end
for i=1:k-1
    plot([chg_pt(i) chg_pt(i)],[P_f_cell{i,1}(end,1)/P_f_MCS P_f_cell{i+1,1}(1,1)/P_f_MCS],'LineWidth',1.5,'Color','m')
end
xlabel('The number of function calls N_{call}','fontsize',11)
ylabel('Normalized P_{f}','fontsize',11)
ylim([0 4])
legend('Estimation','"Exact" solution')
hold off