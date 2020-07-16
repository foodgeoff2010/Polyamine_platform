% cd ../../../GitHub/cobratoolbox;
% initCobraToolbox;

%% Model modification
model = readCbModel('yeast_7.6_COBRA.xml');
% change obj to spermidine production
model = changeObjective(model,'r_2051',1);
% add secretion of 5'-methylthioadenosine
model = addReaction(model,'Sec_s_0303','reactionFormula','s_0303[c_03] -> ');

%% Calculate yield for different pathways

% Pathway1
% add arginine decarboxylase(EC 4.1.1.19) and Agmatinase(EC 3.5.3.11)
model_1 = model;
model_1 = changeRxnBounds(model_1,'r_0817',0,'b');
model_1 = addReaction(model_1,'M1_1','reactionFormula','s_0965[c_03] + s_0794[c_03] -> s_0456[c_03] + agmatine');
model_1 = addReaction(model_1,'M1_2','reactionFormula','agmatine + s_0803[c_03] -> s_1389[c_03] + s_1552[c_03]');
printRxnFormula(model_1,'rxnAbbrList',model_1.rxns(end-1:end),'metNameFlag',1);
sol_1 = optimizeCbModel(model_1,'max','one');

% Pathway2
% add arginine decarboxylase(EC 4.1.1.19), agmatine deaminase(EC 3.5.3.12)
% and N-carbamoylputrescine amidohydrolase(EC 3.5.1.53)
model_2 = model;
model_2 = changeRxnBounds(model_2,'r_0817',0,'b');
model_2 = addReaction(model_2,'M2_1','reactionFormula','s_0965[c_03] + s_0794[c_03] -> s_0456[c_03] + agmatine');
model_2 = addReaction(model_2,'M2_2','reactionFormula','agmatine + s_0803[c_03] -> s_0419[c_03] + N_carbamoylputrescine');
model_2 = addReaction(model_2,'M2_3','reactionFormula','N_carbamoylputrescine + s_0803[c_03] + 2 s_0794[c_03] -> s_1389[c_03] + s_0456[c_03] + s_0419[c_03]');
printRxnFormula(model_2,'rxnAbbrList',model_2.rxns(end-2:end),'metNameFlag',1);
sol_2 = optimizeCbModel(model_2,'max','one');

% Pathway3
% native
model_3 = model;
sol_3 = optimizeCbModel(model_3,'max','one');

% Pathway4
% add carboxyspermidine dehydrogenase (EC 1.5.1.43) and carboxyspermidine
% decarboxylase(EC 4.1.1.-).
model_4 = model;
model_4 = changeRxnBounds(model_4,'r_1001',0,'b');
model_4 = addReaction(model_4,'M4_1','reactionFormula','s_1389[c_03] + s_0978[c_03] + s_0794[c_03] + s_1212[c_03] -> Carboxyspermidine + s_1207[c_03] + s_0803[c_03]');
model_4 = addReaction(model_4,'M4_2','reactionFormula','Carboxyspermidine + s_0794[c_03] -> s_1439[c_03] + s_0456[c_03]');
printRxnFormula(model_4,'rxnAbbrList',model_4.rxns(end-1:end),'metNameFlag',1);
sol_4 = optimizeCbModel(model_4,'max','one');

% Pathway5
% add arginine decarboxylase(EC 4.1.1.19), Agmatinase(EC 3.5.3.11),
% carboxyspermidine dehydrogenase (EC 1.5.1.43) and carboxyspermidine
% decarboxylase(EC 4.1.1.-).
model_5 = model;
model_5 = changeRxnBounds(model_5,'r_0817',0,'b');
model_5 = changeRxnBounds(model_5,'r_1001',0,'b');
model_5 = addReaction(model_5,'M5_1','reactionFormula','s_0965[c_03] + s_0794[c_03] -> s_0456[c_03] + agmatine');
model_5 = addReaction(model_5,'M5_2','reactionFormula','agmatine + s_0803[c_03] -> s_1389[c_03] + s_1552[c_03]');
model_5 = addReaction(model_5,'M5_3','reactionFormula','s_1389[c_03] + s_0978[c_03] + s_0794[c_03] + s_1212[c_03] -> Carboxyspermidine + s_1207[c_03] + s_0803[c_03]');
model_5 = addReaction(model_5,'M5_4','reactionFormula','Carboxyspermidine + s_0794[c_03] -> s_1439[c_03] + s_0456[c_03]');
printRxnFormula(model_5,'rxnAbbrList',model_5.rxns(end-3:end),'metNameFlag',1);
sol_5 = optimizeCbModel(model_5,'max','one');

% Pathway6
% add arginine decarboxylase(EC 4.1.1.19), agmatine deaminase(EC 3.5.3.12),
% N-carbamoylputrescine amidohydrolase(EC 3.5.1.53), carboxyspermidine
% dehydrogenase(EC 1.5.1.43), carboxyspermidine decarboxylase(EC 4.1.1.-).
model_6 = model;
model_6 = changeRxnBounds(model_6,'r_0817',0,'b');
model_6 = changeRxnBounds(model_6,'r_1001',0,'b');
model_6 = addReaction(model_6,'M6_1','reactionFormula','s_0965[c_03] + s_0794[c_03] -> s_0456[c_03] + agmatine');
model_6 = addReaction(model_6,'M6_2','reactionFormula','agmatine + s_0803[c_03] -> s_0419[c_03] + N_carbamoylputrescine');
model_6 = addReaction(model_6,'M6_3','reactionFormula','N_carbamoylputrescine + s_0803[c_03] + 2 s_0794[c_03] -> s_1389[c_03] + s_0456[c_03] + s_0419[c_03]');
model_6 = addReaction(model_6,'M6_4','reactionFormula','s_1389[c_03] + s_0978[c_03] + s_0794[c_03] + s_1212[c_03] -> Carboxyspermidine + s_1207[c_03] + s_0803[c_03]');
model_6 = addReaction(model_6,'M6_5','reactionFormula','Carboxyspermidine + s_0794[c_03] -> s_1439[c_03] + s_0456[c_03]');
printRxnFormula(model_6,'rxnAbbrList',model_6.rxns(end-4:end),'metNameFlag',1);
sol_6 = optimizeCbModel(model_6,'max','one');

% Pathway7
% add arginine decarboxylase(EC 4.1.1.19), agmatine aminopropyl
% transferase(EC 2.5.1.104) and N1-aminopropylagmatine ureohydrolase.
model_7 = model;
model_7 = changeRxnBounds(model_7,'r_0817',0,'b');
model_7 = changeRxnBounds(model_7,'r_1001',0,'b');
model_7 = addReaction(model_7,'M7_1','reactionFormula','s_0965[c_03] + s_0794[c_03] -> s_0456[c_03] + agmatine');
model_7 = addReaction(model_7,'M7_2','reactionFormula','agmatine + s_1420[c_03] -> N1_3_aminopropyl_agmatine + s_0303[c_03] + s_0794[c_03]');
model_7 = addReaction(model_7,'M7_3','reactionFormula','N1_3_aminopropyl_agmatine + s_0803[c_03] -> s_1439[c_03] + s_1552[c_03]');
printRxnFormula(model_7,'rxnAbbrList',model_7.rxns(end-2:end),'metNameFlag',1);
sol_7 = optimizeCbModel(model_7,'max','one');

clear ans;

%% Active reactions and metabolites

[num_rxns_1, num_mets_1] = countActive(model_1,sol_1);
[num_rxns_2, num_mets_2] = countActive(model_2,sol_2);
[num_rxns_3, num_mets_3] = countActive(model_3,sol_3);
[num_rxns_4, num_mets_4] = countActive(model_4,sol_4);
[num_rxns_5, num_mets_5] = countActive(model_5,sol_5);
[num_rxns_6, num_mets_6] = countActive(model_6,sol_6);
[num_rxns_7, num_mets_7] = countActive(model_7,sol_7);

table = [sol_1.f, num_rxns_1, num_mets_1
         sol_2.f, num_rxns_2, num_mets_2
         sol_3.f, num_rxns_3, num_mets_3
         sol_4.f, num_rxns_4, num_mets_4
         sol_5.f, num_rxns_5, num_mets_5
         sol_6.f, num_rxns_6, num_mets_6
         sol_7.f, num_rxns_7, num_mets_7];

clear num_rxns_1 num_mets_1 num_rxns_2 num_mets_2;
clear num_rxns_3 num_mets_3 num_rxns_4 num_mets_4;
clear num_rxns_5 num_mets_5 num_rxns_6 num_mets_6;
clear num_rxns_7 num_mets_7;

clear model;