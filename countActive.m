% Count active reactions and metabolites
function [num_rxns, num_mets] = countActive(model,sol)

flux = sol.x;
idx = flux ~= 0;
activeRxns = model.rxns(idx);
num_rxns = length(activeRxns);

activeMets_tmp = {};
for i = 1:num_rxns
    met_list = model.mets(model.S(:,i) ~= 0);
    
    for j = 1:length(met_list)
        met_id = met_list{j};
        if contains(met_id,'[c')
            met_id_new = met_id(1:6);
        else
            met_id_new = met_id;
        end
        activeMets_tmp = [activeMets_tmp;met_id_new];
    end
end
activeMets = unique(activeMets_tmp);

num_mets = length(activeMets);

