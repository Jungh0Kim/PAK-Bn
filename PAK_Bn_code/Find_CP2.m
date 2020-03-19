% Find critical point

function [CP_hat,find_1_id,find_2_id] = Find_CP2(mu,x_doe,G_doe,N_iniDOE,index_neval,index_eval,Bn_samp_chg)

I_1 = (mu < 0);
I_2 = (G_doe(N_iniDOE+1:end) < 0);
[find_1_id,~] = find(I_1==1);
[find_2_id,~] = find(I_2==1);
if (~isempty(find_1_id) || ~isempty(find_2_id))
    [CP_1,CP_1_id] = max(index_neval(find_1_id));
    [CP_2,CP_2_id] = max(index_eval(find_2_id));
    [~,CP_id] = max([CP_1 CP_2]);
    if CP_id==1
        CP_hat(1,:) = Bn_samp_chg(find_1_id(CP_1_id),:);
    elseif CP_id==2
        CP_hat(1,:) = x_doe(N_iniDOE + find_2_id(CP_2_id),:);
    end
end

end % function end


