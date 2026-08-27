% 2026.05.13 - Author: Lucas Santana Souza 
% 
% Aim: save data for max growth in different conditions: 
%      [ancestor of host  alone in sole env]; 
%      [ancestor of guest alone in sole env]; 
%      [pair with equal growth], 
%      [ancestors of host  in shared env]
%      [ancestors of guest in shared env]
%
%% Structure of 'run_AloneNonShared_cfe(ehmodel1 ,ne ,nic ,nio)':
% [ S_ext ]  [>] [(Host's rhs_ext_lb) ]-> REGION b1
% [ S_ext ]  [<] [(Host's rhs_ext_ub) ]-> REGION b2  
% [ S_intC ] [=] [0]                   -> REGION b3 
% [ S_intO ] [=] [0]                   -> REGION b4 
%
%% Structure of 'create_AloneShared_cfe (ehmodel1 ,ehmodel2 ,ne ,nic ,nio)':
% [ S_ext ]   [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [ S_ext ]   [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [ S_intC ]  [=] [0]                                           -> REGION b3 
% [ S_intO ]  [=] [0]                                           -> REGION b4 
%
%% Structure of 'create_holobiont_default_rhs_cfe (ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu)':
%   HOST          GUEST
% [  S_ext ] | [      0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [  S_ext ] | [      0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [    0   ] | [ S_unmapped ] [=] [0]                                           -> REGION b3  
% [ S_intC ] | [ S_ext2intC ] [=] [0]                                           -> REGION b4 
% [    0   ] | [   S_intC   ] [=] [0]                                           -> REGION b5  
% [    0   ] | [   S_intO   ] [=] [0]                                           -> REGION b6  
% [ S_intO ] | [     0      ] [=] [0]                                           -> REGION b7  
% [    1   ] | [    -1      ] [=] [0]   -> REGION b8 for s.t (Host's growth) = (Guest's growth)
%
%
% 'S_ext' -> Host's compartment only contain external metabolites that can be mapped
% 'S_intC' -> Contain intracellular metabolites from Citoplasm
% 'S_intO' -> Contain intracellular metabolites from Other Compartments:
%             (e.g. Peroxisome, Mitochondrion, Endoplasmic reticulum, Lipid particle, Nucleus, and Golgi)
% 'S_ext2intC' ->  guest's ext compartment contain external metabolites that can be mapped    to host's intC
% 'S_unmapped' ->  guest's ext compartment contain external metabolites that cannot be mapped to host's intC 
%
% PS.: 
% To avoid transparency violation in parallel runs, I replaced 'clear' with '= []'

clear all
clc
format compact

%% Running criteria

tic;

%importPairsOrGenerateHere = 'importPairs'
importPairsOrGenerateHere = 'GenerateHere' 

%runningIn = 'hpc2n'
runningIn = 'desktop'

% Saving data Sufix
sufix_database = '_CarveFungiEndo'; 

rng(0) %initializes Mersenne Twister algorithm with seed 0

sampleSizeUsed= 2000000;


%% Gurobi settings  
addpath 'C:\gurobi1001\win64\matlab'
params = struct();
params.OutputFlag = 0;
params.FeasibilityTol= 1e-9;
params.OptimalityTol = 1e-9; % Default value in gurobi is: 1e-6

%% Defining upload directory
% Define the username                                                      % change dependening on computer
%username = 'lsant';
username =  'lucas';
%username = 'lusa4312';

% Define where data is uploaded from: Onedrive (cloud) or local  
%cloud_local = '/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
%cloud_local = '/Desktop/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL

% Define the data used
dataUsed = '/CarveFungiEndo_rb'

% Define the input directory path for INPUT: metabolic model
%filedir='C:\Users\lusa4312\Documents\ProkaryoteEndosymbiosis-main\renamingAGORA'; % alternative code for the one bellow 
filedirGEM = ['C:/Users/'  username  cloud_local  dataUsed ]


%% Get a list of all files in the directory
files = dir(fullfile(filedirGEM, '*.mat'));

% Count the number of ".mat" files
num_mat_files = length(files)



%% Defyning pairs to use

% Define where data is uploaded from: Onedrive (cloud) or local  
if strcmp (username,'lucas' )   
    cloud_local = '/Desktop/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
end
if strcmp (username, 'lusa4312' )   
    cloud_local = '/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
end



% Count the number of ".mat" files
if strcmp (importPairsOrGenerateHere, 'importPairs')
    
    % Define the data used
    sufixdir = '/growthResults_cfe/data_in_format_to_analyse';
    
    % Define the input directory path for INPUT: idH, idG 
    dir_UploadIDS = ['C:/Users/'  username  cloud_local  sufixdir ]
    
    % Uploading IDs
    input_file_name = ['IdRow_IdH_IdG_cfe_Sampled'  num2str(sampleSizeUsed) '.csv'];            
    hId_gId = fullfile(dir_UploadIDS, input_file_name);   
    rowID_hID_gID_sampled = readmatrix(hId_gId);       % Data structure: 1st-column is hostID, 2nd-column is guestID  
    
    simulationSize =  size(rowID_hID_gID_sampled, 1)
end

if strcmp (importPairsOrGenerateHere, 'GenerateHere')    % Running this section takes about 4min
  % Defining aspects used to create all possible combination of pairs

  position  = 1;
  n_rows      = num_mat_files*num_mat_files % Total # rows for all combination of pairs
  hID_gID_all = zeros (n_rows,2); % Pre-empting array that will contain all pair combinations possible

  % Creating all possible combination of host-guest pairs
  for col1 = 1:num_mat_files
      for col2 = 1:num_mat_files
          hID_gID_all (position,1) = col1;
          hID_gID_all (position,2) = col2;
          position = position + 1;
      end
  end
  position =1;

  rowID = [1:(num_mat_files*num_mat_files)]';

  % Combining the data
  rowID_hID_gID = [rowID    hID_gID_all];

    %% Upload ID of sp without growth
    
    % Define the input directory path for INPUT: spId cannot grow alone in NonShared Env  
    if strcmp (username,'lucas' )   
        dir_UploadIDS = ['C:\Users\lucas\Desktop\OneDrive\ProkaryoteEndosymbiosis-main\growthResults_cfe\data_in_format_to_analyse' ]
    end
    if strcmp (username, 'lusa4312' )   
        dir_UploadIDS = ['C:\Users\lusa4312\OneDrive\ProkaryoteEndosymbiosis-main\growthResults_cfe\data_in_format_to_analyse' ]
    end
    
    
    % Uploading IDs
    input_file_name = ['v_SpID_wAloneNonSharedEnv_NoGrowth_cfe_default_rhs_CarveFungiEndo.csv']            
    SpId_noGrowth = fullfile(dir_UploadIDS, input_file_name);   
    remove_values = readmatrix(SpId_noGrowth);        
    
    col = 2;
    rowID_hID_gID(ismember(rowID_hID_gID(:,col), remove_values), :) = [];
    
    col = 3;
    rowID_hID_gID(ismember(rowID_hID_gID(:,col), remove_values), :) = [];

    
    % Sampling IDs that will be used in simulation
    
    % Defining size of the sample, it must be: simulationSize < (number of GEMs)^2
    simulationSize = sampleSizeUsed                                               % HERE critical place to change
    n_rows = size (rowID_hID_gID,1) % This is give the number of all possible pairs when excluding all sp that cannot grow alone in default NonShared env
    
    sampled_IDs = randsample(n_rows ,simulationSize); %returns values sampled uniformly at random, without replacement, from the intergers 1 to n_rows.
    sampled_IDs = sort(sampled_IDs);
    
    rowID_hID_gID_sampled = rowID_hID_gID(sampled_IDs,:); %list of host-guest pairs that will be used
    hID_gID_sampled =  rowID_hID_gID_sampled (:,[2 3]);
    


    % Saving the IDS
        
    % Define the output directory path
    if strcmp (runningIn, 'desktop')
        SavingPathName = ['C:\Users\' username cloud_local '\growthResults_cfe\growthResultsTestID']; %
    elseif strcmp (runningIn, 'hpc2n')   
        SavingPathName = '/pfs/proj/nobackup/fs/projnb10/hpc2n2023-112/lusa4312/Documents/ProkaryoteEndosymbiosis-main/growthResults_cfe/growthResultsTestID';
    end
    
    
    % Prefix used in saved name:
    hID_gID_prefix       =       ['IdH_IdG_cfe_Sampled'  num2str(simulationSize) ];
    rowID_hID_gID_prefix = ['IdRow_IdH_IdG_cfe_Sampled'  num2str(simulationSize) ];
    
    %Save sampled IDs
    writematrix(       hID_gID_sampled  ,fullfile(SavingPathName ,      [hID_gID_prefix  sufix_database '.csv']));
    writematrix( rowID_hID_gID_sampled  ,fullfile(SavingPathName ,[rowID_hID_gID_prefix  sufix_database '.csv']));


end




%% Separating sampled values to use in loop 
rowIds = rowID_hID_gID_sampled (:,1);    
hIds   = rowID_hID_gID_sampled (:,2);
gIds   = rowID_hID_gID_sampled (:,3);
 
%% Setting arrays that storage the max growth rate
gv_ancestral_alone_nonSharedEnv_h = zeros (simulationSize ,1);
gv_ancestral_alone_nonSharedEnv_g = zeros (simulationSize ,1);
gv_pair_st_equal_growth           = zeros (simulationSize ,1);
gv_ancestral_alone_SharedEnv_h    = zeros (simulationSize ,1);
gv_ancestral_alone_SharedEnv_g    = zeros (simulationSize ,1);

gv_h_name = strings(simulationSize,1);
gv_g_name = strings(simulationSize,1);





%% Finding parameters to use within loop
% Step1: Replaced EVAL with LOAD to load the host model, avoiding "transparency violation error"
host_model_path = fullfile(filedirGEM, ['model1.mat']);
ehmodel1 = load(host_model_path, 'metabolic_model');
ehmodel1 = ehmodel1.metabolic_model;


% Step2: Parameters used in function to create the GEM structure 
ne  = size(ehmodel1.S_ext      ,1); % ne  -> # of extracelular mapped metabolites
nic = size(ehmodel1.S_intC     ,1); % nic -> # of metabolites in citoplasm
nio = size(ehmodel1.S_intO     ,1); % nio -> # of metabolites in other compartments
nu  = size(ehmodel1.S_unmapped ,1); % nu  -> # of extracelular unmapped metabolites

% Step3: Cleaning
ehmodel1 =[];


%% Initialize a parallel pool
parpool;

% Calculation
parfor i = 1:simulationSize
    
    % Specify which sample will be used 
    %rowNum = rowIds(i); % If I attempt to index 'rowNum' instead of 'i', the parfor does not run.
    ind1   = hIds (i); % integer for a metabolic model for host
    ind2   = gIds (i); % integer for a metabolic model for host

    %% compute growth rates for ancestrals in own environs
    
    % Load the host model
    % Change 1: Replaced EVAL with LOAD to load the host model, avoiding "transparency violation error"
    h_model_path = fullfile(filedirGEM, ['model', num2str(ind1), '.mat']);
    ehmodel1 = load(h_model_path, 'metabolic_model');
    ehmodel1 = ehmodel1.metabolic_model;
    
    % Outputs host's ancestral growth rate s.t. growing alone in its own env (non-shared env)
    [resulth ,hmodel] = run_AloneNonShared_cfe(ehmodel1 ,ne ,nic ,nio);
            
    % Stores the max growth rates of ancestors growing alone in nonShared env
    gv_ancestral_alone_nonSharedEnv_h(i)= abs(resulth.objval); 

    % Sp name
    gv_h_name(i) = string(ehmodel1.description);


    % Load the guest (endo) model
    % Change 2: Replaced EVAL with LOAD to load the endo model, avoiding "transparency violation error"
    g_model_path = fullfile(filedirGEM, ['model', num2str(ind2), '.mat']);
    ehmodel2 = load(g_model_path, 'metabolic_model');
    ehmodel2 = ehmodel2.metabolic_model;

    % Outputs guest's ancestral growth rate s.t. growing alone in its own env (non-shared env)
    [resultg, gmodel] = run_AloneNonShared_cfe(ehmodel2 ,ne ,nic ,nio)
    
    % Stores the max growth rates of ancestors growing alone in nonShared env
    gv_ancestral_alone_nonSharedEnv_g(i)= abs(resultg.objval); 

    % Sp name
    gv_g_name(i) = string(ehmodel2.description);

    %% Create holobiont model s.t. (host's growth rate) = (guest's growth rate) 
    %  pairmodel = create_holobiont_default_rhs_Ng (ehmodel1,ehmodel2 ,ne ,ni, nu ,Ng);   
    pairmodel = create_holobiont_default_rhs_cfe (ehmodel1,ehmodel2,ne ,nic ,nio ,nu);   
    
    % Outputs host's growth rate s.t. (host's growth rate) = (guest's growth rate)
    resultpairhe = gurobi(pairmodel,params);

    % Check if feasible & extract growth rate
    if ~strcmp(resultpairhe.status,'OPTIMAL')
        resultpairhe=0;
    else
        resultpairhe=abs(resultpairhe.objval);
    end

    %% Grow alone in Shared Environment
    
    %%  Maximize growth of alone in shared env. (same ID of host)       
    modelHAS =  create_AloneShared_cfe (ehmodel1 ,ehmodel2 ,ne ,nic ,nio);
    
    % Outputs alone in shared environment (used same ID of host)
    gurobi_result_HAS = gurobi(modelHAS ,params);
    wHAS  = [];
    % check if feasible
    if ~strcmp(gurobi_result_HAS.status ,'OPTIMAL')
        wHAS = 0;                           
    else
        wHAS = abs(gurobi_result_HAS.objval);
    end
    
    %%  Maximize growth of alone in shared env. (same ID of guest)       
    modelGAS =  create_AloneShared_cfe (ehmodel2 ,ehmodel1 ,ne ,nic ,nio);
    
    % Outputs alone in shared environment (used same ID of host)
    gurobi_result_GAS = gurobi(modelGAS ,params);
    wGAS  = [];
    % check if feasible
    if ~strcmp(gurobi_result_GAS.status ,'OPTIMAL')
        wGAS = 0;                           
    else
        wGAS = abs(gurobi_result_GAS.objval);
    end

    %% Stores the max growth rates
    gv_pair_st_equal_growth (i) = resultpairhe;
    gv_ancestral_alone_SharedEnv_h (i) = wHAS;
    gv_ancestral_alone_SharedEnv_g (i) = wGAS;

    %% Clean structs
    ehmodel2           = [];
    pairmodel          = [];
    resultpairhe       = [];
    modelHAS           = [];
    modelGAS           = [];
    wHAS               = [];
    wGAS               = [];
    gurobi_result_GAS  = []; 
    gurobi_result_HAS  = [];

    %end   % internal for loop 

    %% Clean structs
    metabolic_model = [];
    ehmodel1        = [];
    resultg         = [];
    hmodel          = [];
    gmodel          = [];

end

% Close the parallel pool
delete(gcp('nocreate'));

%% Defining download directory

% Define where data is uploaded from: Onedrive (cloud) or local  
%cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL
if strcmp (username, 'lucas')   
    cloud_local = '/Desktop/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
end
if strcmp (username, 'lusa4312')   
    cloud_local = '/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
end

% Define the output directory path
%SavingPathName = 'C:\Users\lusa4312\OneDrive\ProkaryoteEndosymbiosis-main\growthResults_cfe\growthResultsTestP'; % alternative code for the one bellow 
SavingPathName = ['C:\Users\' username cloud_local '\growthResults_cfe\growthResultsTestP']; % P:= growth in sync

%% Saving data

% Prefix used in saved name:
alone_h_NonShared_prefix            = 'v_rowID_hID_gID_wAloneHNonSharedEnv_hName_cfe_default_rhs';
alone_g_NonShared_prefix            = 'v_rowID_hID_gID_wAloneGNonSharedEnv_gName_cfe_default_rhs';
gv_pair_st_equal_growth_prefix      = 'v_rowID_hID_gID_wSync_cfe_default_rhs' ;
ancestral_alone_h_SharedEnv_prefix  = 'v_rowID_hID_gID_wAloneHSharedEnv_cfe_default_rhs';
ancestral_alone_g_SharedEnv_prefix  = 'v_rowID_hID_gID_wAloneGSharedEnv_cfe_default_rhs';

% Saving data Sufix
sufix_database = '_CarveFungiEndo';  

% including rowId, hID and gID to the file saved
rowID_hID_gID_alone_nonSharedEnv_h        = [rowID_hID_gID_sampled   gv_ancestral_alone_nonSharedEnv_h   gv_h_name];
rowID_hID_gID_alone_nonSharedEnv_g        = [rowID_hID_gID_sampled   gv_ancestral_alone_nonSharedEnv_g   gv_g_name];
rowID_hID_gID_sync                        = [rowID_hID_gID_sampled   gv_pair_st_equal_growth];
rowID_hID_gID_ancestral_alone_SharedEnv_h = [rowID_hID_gID_sampled   gv_ancestral_alone_SharedEnv_h];
rowID_hID_gID_ancestral_alone_SharedEnv_g = [rowID_hID_gID_sampled   gv_ancestral_alone_SharedEnv_g];

% gv_h_name = ehmodel1.description
% gv_g_name = ehmodel2.description


% Saving files
writematrix(rowID_hID_gID_alone_nonSharedEnv_h        ,fullfile(SavingPathName ,[alone_h_NonShared_prefix           '_Sampled'  num2str(simulationSize)  sufix_database '.csv' ]));
writematrix(rowID_hID_gID_alone_nonSharedEnv_g        ,fullfile(SavingPathName ,[alone_g_NonShared_prefix           '_Sampled'  num2str(simulationSize)  sufix_database '.csv' ]));
writematrix(rowID_hID_gID_sync                        ,fullfile(SavingPathName ,[gv_pair_st_equal_growth_prefix     '_Sampled'  num2str(simulationSize)  sufix_database '.csv' ]));
writematrix(rowID_hID_gID_ancestral_alone_SharedEnv_h ,fullfile(SavingPathName ,[ancestral_alone_h_SharedEnv_prefix '_Sampled'  num2str(simulationSize)  sufix_database '.csv' ]));
writematrix(rowID_hID_gID_ancestral_alone_SharedEnv_g ,fullfile(SavingPathName ,[ancestral_alone_g_SharedEnv_prefix '_Sampled'  num2str(simulationSize)  sufix_database '.csv' ]));



toc



