% Find critical point

function [CP_hat,find_1_id] = Find_CP1(mu,index_neval,Bn_samp)

I_1 = (mu < 0);
[find_1_id,~] = find(I_1==1);
if ~isempty(find_1_id)
    [~,CP_id] = max(index_neval(find_1_id));
    CP_hat(1,:) = Bn_samp(find_1_id(CP_id),:);
end

end % function end


