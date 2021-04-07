% cd ../../../GitHub/cobratoolbox;
% initCobraToolbox;

%% Model modification
model = readCbModel('yeast_7.6_COBRA.xml');
% change obj to spermidine production
model = changeObjective(model,'r_2051',1);
% add secretion of 5'-methylthioadenosine
model = addReaction(model,'Sec_s_0303','reactionFormula','s_0303[c_03] -> ');

% printRxnFormula(model,'rxnAbbrList',model.rxns);

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

clear model;

%% New tasks
model = readCbModel('yeast_7.6_COBRA.xml');

% Pathway8
% Theoretical yield of Spermine
model_8 = model;
% change obj to spermine exchange
model_8 = changeObjective(model_8,'r_2052',1);
sol_8 = optimizeCbModel(model_8,'max','one');

% Pathway9
% Theoretical yield of Thermospermine
% add thermospermine synthase(EC 2.5.1.79) and export of thermospermine
model_9 = model;
model_9 = addReaction(model_9,'M9_1','reactionFormula','s_1420[c_03] + s_1439[c_03] -> s_0303[c_03] + thermospermine + s_0794[c_03]');
model_9 = addReaction(model_9,'M9_2','reactionFormula','thermospermine + s_0796[c_06] -> s_0794[c_03] + thermospermine_e');
model_9 = addReaction(model_9,'M9_3','reactionFormula','thermospermine_e -> ');
printRxnFormula(model_9,'rxnAbbrList',model_9.rxns(end-2:end),'metNameFlag',1);
model_9 = changeObjective(model_9,'M9_3',1);
sol_9 = optimizeCbModel(model_9,'max','one');

% Pathway10
% Theoretical yield of sym-homospermidine (biosynthesis I)
% add homospermidine synthase(EC 2.5.1.45) and export of sym-homospermidine
% and trimethylenediamine
model_10 = model;
model_10 = addReaction(model_10,'M10_1','reactionFormula','s_1439[c_03] + s_1389[c_03] -> sym-homospermidine + s_1526[c_03]');
model_10 = addReaction(model_10,'M10_2','reactionFormula','sym-homospermidine -> sym-homospermidine_e');
model_10 = addReaction(model_10,'M10_3','reactionFormula','sym-homospermidine + s_0796[c_06] -> s_0794[c_03] + sym-homospermidine_e');
model_10 = addReaction(model_10,'M10_4','reactionFormula','sym-homospermidine_e -> ');
model_10 = addReaction(model_10,'M10_5','reactionFormula','s_1526[c_03] -> trimethylenediamine_e');
model_10 = addReaction(model_10,'M10_6','reactionFormula','s_1526[c_03] + s_0796[c_06] -> s_0794[c_03] + trimethylenediamine_e');
model_10 = addReaction(model_10,'M10_7','reactionFormula','trimethylenediamine_e -> ');
printRxnFormula(model_10,'rxnAbbrList',model_10.rxns(end-6:end),'metNameFlag',1);
model_10 = changeObjective(model_10,'M10_4',1);
sol_10 = optimizeCbModel(model_10,'max','one');

% Pathway11
% Theoretical yield of sym-homospermidine (biosynthesis II)
% add homospermidine synthase(EC 2.5.1.44) and export of sym-homospermidine
model_11 = model;
model_11 = addReaction(model_11,'M11_1','reactionFormula','2 s_1389[c_03] -> sym-homospermidine + s_0419[c_03]');
model_11 = addReaction(model_11,'M11_2','reactionFormula','sym-homospermidine -> sym-homospermidine_e');
model_11 = addReaction(model_11,'M11_3','reactionFormula','sym-homospermidine + s_0796[c_06] -> s_0794[c_03] + sym-homospermidine_e');
model_11 = addReaction(model_11,'M11_4','reactionFormula','sym-homospermidine_e -> ');
printRxnFormula(model_11,'rxnAbbrList',model_11.rxns(end-3:end),'metNameFlag',1);
model_11 = changeObjective(model_11,'M11_4',1);
sol_11 = optimizeCbModel(model_11,'max','one');

clear ans;

%% Active reactions and metabolites

[num_rxns_1, num_mets_1] = countActive(model_1,sol_1);
[num_rxns_2, num_mets_2] = countActive(model_2,sol_2);
[num_rxns_3, num_mets_3] = countActive(model_3,sol_3);
[num_rxns_4, num_mets_4] = countActive(model_4,sol_4);
[num_rxns_5, num_mets_5] = countActive(model_5,sol_5);
[num_rxns_6, num_mets_6] = countActive(model_6,sol_6);
[num_rxns_7, num_mets_7] = countActive(model_7,sol_7);
[num_rxns_8, num_mets_8] = countActive(model_8,sol_8);
[num_rxns_9, num_mets_9] = countActive(model_9,sol_9);
[num_rxns_10, num_mets_10] = countActive(model_10,sol_10);
[num_rxns_11, num_mets_11] = countActive(model_11,sol_11);

table = [sol_1.f, num_rxns_1, num_mets_1
         sol_2.f, num_rxns_2, num_mets_2
         sol_3.f, num_rxns_3, num_mets_3
         sol_4.f, num_rxns_4, num_mets_4
         sol_5.f, num_rxns_5, num_mets_5
         sol_6.f, num_rxns_6, num_mets_6
         sol_7.f, num_rxns_7, num_mets_7
         sol_8.f, num_rxns_8, num_mets_8
         sol_9.f, num_rxns_9, num_mets_9
         sol_10.f, num_rxns_10, num_mets_10
         sol_11.f, num_rxns_11, num_mets_11];

clear num_rxns_1 num_mets_1 num_rxns_2 num_mets_2;
clear num_rxns_3 num_mets_3 num_rxns_4 num_mets_4;
clear num_rxns_5 num_mets_5 num_rxns_6 num_mets_6;
clear num_rxns_7 num_mets_7 num_rxns_8 num_mets_8;
clear num_rxns_9 num_mets_9 num_rxns_10 num_mets_10;
clear num_rxns_11 num_mets_11;

%%
% model_tmp = model_11;
% sol_tmp = sol_11;
% idx = sol_tmp.x ~= 0;
% z1_tmp = model_tmp.rxns(idx);
% z2_tmp = printRxnFormula(model_tmp,'rxnAbbrList',z1_tmp,'metNameFlag',1);
% z3_tmp = sol_tmp.x(idx);
