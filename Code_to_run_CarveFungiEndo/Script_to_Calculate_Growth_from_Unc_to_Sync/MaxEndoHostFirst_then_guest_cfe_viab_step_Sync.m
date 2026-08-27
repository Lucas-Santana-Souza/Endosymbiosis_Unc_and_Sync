% 2026.05.27 - Author: Lucas Santana Souza 
%
% I used within this code:
% create_pair_step1_ofMHF_cfe_viab          (ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,minGrowth)
% create_pair_step2_ofMHF_cfe               (ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,result_pair12_step1_ofMHF_firstM1)
% create_pair_step1_ofMHF_cfe_viab_step_Sync(ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,minGrowth ,result_pair12_step2_ofMHF_firstM1 ,stepSize  ,Sync_Fitness);
% create_pair_step2_ofMHF_cfe               (ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,result_pair12_step1_ofMHF_firstM1); 
%
% Aim: to calculate growth rates for hosts having increasing investment (stepSize) in guest's growth
%      which will be used to calculate trade-off curves. 
%  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START of Initialization at starting time ( before while loop, t=1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% create_pair_step1_ofMHF_cfe_viab (ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,minGrowth); 
% [growth rate of the Host maximized first, s.t. (guest's growth) > minGrowth]
%   HOST          GUEST
% [  S_ext ] | [      0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [  S_ext ] | [      0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [    0   ] | [ S_unmapped ] [=] [0]                                           -> REGION b3  
% [ S_intC ] | [ S_ext2intC ] [=] [0]                                           -> REGION b4 
% [    0   ] | [   S_intC   ] [=] [0]                                           -> REGION b5  
% [    0   ] | [   S_intO   ] [=] [0]                                           -> REGION b6  
% [ S_intO ] | [      0     ] [=] [0]                                           -> REGION b7  
% [    0   ] | [      1     ] [>] [minGrowth]   -> REGION b8 for s.t (guest's growth) > minGrowth
%
%% create_pair_step2_ofMHF_cfe (ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,result_pair12_step1_ofMHF_firstM1); 
% [growth rate of the Guest maximized second, s.t. (host max first)]
%   HOST          GUEST
% [  S_ext ] | [      0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [  S_ext ] | [      0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [    0   ] | [ S_unmapped ] [=] [0]                                           -> REGION b3  
% [ S_intC ] | [ S_ext2intC ] [=] [0]                                           -> REGION b4 
% [    0   ] | [   S_intC   ] [=] [0]                                           -> REGION b5  
% [    0   ] | [   S_intO   ] [=] [0]                                           -> REGION b6  
% [ S_intO ] | [      0     ] [=] [0]                                           -> REGION b7  
% [    1   ] | [      0     ] [>] [result_pairij_step1_ofMHF_firstMi - ErrorTolerance] -> REGION b8 for constrain in the growth rate when guest is maximized, s.t. host max first
%
%
% OBS: The solution of 'create_pair_step2_ofMHF_uei_Ng'  ( which is  result_pair12_step2_ofMHF_firstM1) will then be used within the while loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END of  of initialization at starting time ( before while loop)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START of  Initialization of while loop time ( t>1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% create_pair_step1_ofMHF_cfe_viab_step_Sync(ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,minGrowth ,result_pair12_step2_ofMHF_firstM1 ,stepSize  ,Sync_Fitness);
% [growth rate of the Host maximized first, s.t. (guest's growth) > minGrowth]
%    HOST           GUEST
% [  S_ext ] | [      0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )] -> REGION b1
% [  S_ext ] | [      0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )] -> REGION b2  
% [    0   ] | [ S_unmapped ] [=] [0]                                            -> REGION b3  
% [ S_intC ] | [ S_ext2intC ] [=] [0]                                            -> REGION b4 
% [    0   ] | [   S_intC   ] [=] [0]                                            -> REGION b5  
% [    0   ] | [   S_intO   ] [=] [0]                                            -> REGION b6  
% [ S_intO ] | [      0     ] [=] [0]                                            -> REGION b7  
% [    0   ] | [      1     ] [>] [minGrowth]                                    -> REGION b8 for s.t (guest's growth) > minGrowth
% [    0   ] | [      1     ] [>] [result_pair12_step2_ofMHF_firstM1 + stepSize] -> REGION b9 for the constrain: (g2Endo's growth) > (result_pair12_step2_ofMHF_firstM1 + stepSize)
% [    1   ] | [      0     ] [>] [Sync_Fitness - ErrorTolerance]                -> REGION b10 for the constrain: h1Endo's growth > Sync_Fitness 

%
%% create_pair_step2_ofMHF_cfe (ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,result_pair12_step1_ofMHF_firstM1); 
% {growth rate of the Guest maximized second, s.t. (host max first) }
%   HOST          GUEST
% [  S_ext ] | [      0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [  S_ext ] | [      0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [    0   ] | [ S_unmapped ] [=] [0]                                           -> REGION b3  
% [ S_intC ] | [ S_ext2intC ] [=] [0]                                           -> REGION b4 
% [    0   ] | [   S_intC   ] [=] [0]                                           -> REGION b5  
% [    0   ] | [   S_intO   ] [=] [0]                                           -> REGION b6  
% [ S_intO ] | [      0     ] [=] [0]                                           -> REGION b7  
% [    1   ] | [      0     ] [>] [result_pairij_step1_ofMHF_firstMi - ErrorTolerance] -> REGION b8 for constrain in the growth rate when guest is maximized, s.t. host max first
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END of  Initialization of while loop time ( t>1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%              ehmodel1 (or hemodel1) -> microbe that is the host
%              ehmodel2 (or hemodel2) -> microbe that is the guest
%                            'S_ext'  -> Host's compartment only contain external metabolites that can be mapped
%                            'S_intC' -> Host's/Guest's compartment only contain internal metabolites from Citoplasm
%                            'S_intO' -> Host's/Guest's compartment only contain internal metabolites from Other compartments (e.g. nucleos, golgi…)
%                        'S_ext2intC' -> metabolictes from guest's extracellular that can    be mapped to host's intC
%                       'S_unmapped'  -> metabolictes from guest's extracellular that cannot be mapped to host's intC 
%                                 ne  -> # of extracelular mapped metabolites
%                                 nic -> # of metabolites within S_intC
%                                 nio -> # of metabolites within S_intO
%                                 nu  -> # of extracelular unmapped metabolites
%                          'stepSize' -> increment for G's growth
% 'result_pair12_step2_ofMHF_firstM1' -> guest max 2nd in endosymbiosis
%

clear all
clc
format compact
%%
disp ('Sim. Started')
t0 = tic;

%% Defining where I am running
runningIn = 'desktop';
%runningIn = 'hpc2n';

%% Defyning if I am running a test simulation (a subset of metabolic models)
%typeOfRun = '3cases'
%typeOfRun = 'ImportCasesReadyToRun'
typeOfRun = 'ImportCases_e_SampleHere'

%% Gurobi settings  
if strcmp (runningIn, 'desktop')
      addpath 'C:\gurobi1001\win64\matlab'                                                                     
elseif strcmp (runningIn, 'hpc2n')   
      %addpath 'C:\gurobi1001\win64\matlab';  %to run in cluster this must be commented;
end

params = struct();
params.OutputFlag = 0;
params.FeasibilityTol=1e-9;
params.OptimalityTol  = 1e-9; % Default value in gurobi is: 1e-6

%%  thres_min_growth is used in the second step of the optimization
thres_min_growth = 1e-8;  %0.001;                                          % this is to not allow the second optimization to occur if the growth rate of the first optimized is below this value

% Constraint applyied within create pair function: specificaly: first step (max first) by making: seconds' growth > minGrowth
minGrowth = 0 %0.001

% To be used in prefix of saved name:
%Extract the value that max 1st is constrained during its optimization. Constraint: (max 2nd > minGrowth) 
strNumber = num2str(minGrowth);  % Convert number to string
% Find the position of the decimal point
dotIndex = strfind(strNumber, '.');
% Extract everything after the decimal point
minGrowthStr = strNumber(dotIndex+1:end);


if length(minGrowthStr) == 0 
   minGrowthStr = strNumber;  
end

%% Define the username                                                      % change dependening on computer
username =  'lucas'
%username = 'lusa4312'

% Define the collection used
collection = 'CarveFungiEndo_rb'

rng(0);  % Setting seed for reproducibility (Use any integer value to set a random seed)


%% DEFINE Database for upload growth rates already calculated
% Define where data is uploaded from: Onedrive (cloud) or local  
%cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL

if strcmp (username, 'lusa4312')
    cloud_local = '/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
end

if strcmp (username, 'lucas')
    cloud_local = '/Desktop/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
end

database = '/growthResults_cfe/';  % Define the location of input: growth of specified collection
sufix_database = '_CarveFungiEndo';                            % Define the sufix of input file







%% Upload hostID and guestID 


% Count the number of ".mat" files
if strcmp (typeOfRun, '3cases')
    hIds   =  [7 ; 8]         
    gIds   =  [2558 ; 1475]                    
    vec_P  =  [ 9.70052267459332; 10.1890590964123 ]   
end


if strcmp (typeOfRun, 'ImportCasesReadyToRun')

    % Define the directory path for INPUT: growth rate and IDs
    %dir_pathGrowth = 'C:\Users\lusa4312\OneDrive\ProkaryoteEndosymbiosis-main\growthResults_cfe\7HDmP_HASmP_1HpGm_VH1G2Endo_Vp_1300';
    dir_pathGrowth  = ['C:/Users/'  username  cloud_local  database '7HDmP_HASmP_1HpGm_VH1G2Endo_Vp_1300'] 
    
    
    % Upload CSV file for pair fitness assuming syncronization.
    input_file_name = ['IdH_IdG_NameH_NameG_WSync_for_7HdmP_HASmP_1HpGm_VH1G2Endo_VP_viab0_cfe_default_rhs_SubSampled1300_CarveFungiEndo.csv']            
    P = fullfile(dir_pathGrowth, input_file_name)
    T_P = readtable(P);

    hIds   = T_P{:, 1};  % 2026.05.25 - using '{}' extracts the raw values from column and converts table format into a array.
    gIds   = T_P{:, 2};  % 2026.05.25 
    vec_P  = T_P{:, 5};  % 2026.05.25 only using the values of Wsync
end


if strcmp (typeOfRun, 'ImportCases_e_SampleHere')
    % Define the directory path for INPUT: growth rate and IDs
    %dir_pathGrowth = 'C:\Users\lusa4312\OneDrive\ProkaryoteEndosymbiosis-main\growthResults_cfe\7HDmP_HASmP_1HpGm_VH1G2Endo_Vp_1300';
    dir_pathGrowth  = ['C:/Users/'  username  cloud_local  database '7HDmP_HASmP_1HpGm_VH1G2Endo_Vp_c'] 
    
    
    % Upload CSV file for pair fitness assuming syncronization.
                  %    'IdH_IdG_NameH_NameG_WSync_for_7HdmP_HASmP_1HpGm_VH1G2Endo_VP_viab0_cfe_default_rhs_Sampled2000000_CarveFungiEndo'
    input_file_name = ['IdH_IdG_NameH_NameG_WSync_for_7HdmP_HASmP_1HpGm_VH1G2Endo_VP_viab0_cfe_default_rhs_Sampled2000000_CarveFungiEndo_c.csv']  
                       % This file contains all the cases in which 7HdmP_HASmP_1HpGm_VH1G2Endo_VP is true out of a random pickup of pairs of size 2 000 000,
                       

    P = fullfile(dir_pathGrowth, input_file_name)
    T_P = readtable(P);

    % Sample rows without replacement
    sample_size = 8%000 % Obs: sample_size < height(T_P)
    idx_sample = randperm(height(T_P), sample_size);
    T_P_sample = T_P(idx_sample, :);

    hIds   = T_P_sample{:, 1};  % 2026.05.25 - using '{}' extracts the raw values from column and converts table format into a array.
    gIds   = T_P_sample{:, 2};  % 2026.05.25 
    vec_P  = T_P_sample{:, 5};  % 2026.05.25 only using the values of Wsync


    if strcmp (username, 'lusa4312')
        cloud_local = '/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
    end
    
    if strcmp (username, 'lucas')
        cloud_local = '/Desktop/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
    end

    SavingPathName = ['C:\Users\' username cloud_local '\growthResults_cfe\growthResultsTestH']
    %SavingPathName = 'C:\Users\lusa4312\OneDrive\ProkaryoteEndosymbiosis-main\growthResults_cfe\growthResultsTestH';

    writetable(T_P_sample ...
          , fullfile(SavingPathName, ['IdH_IdG_NameH_NameG_wSync_for_7HdmP_HASmP_1HpGm_SubSampled' num2str(sample_size)  '_CarveFungiEndo.csv']) ...
          , 'WriteVariableNames', false);

    writematrix(T_P_sample{:, [1 2]} ...
      , fullfile(SavingPathName, ['IdH_IdG_for_7HdmP_HASmP_1HpGm_SubSampled' num2str(sample_size)  '_CarveFungiEndo.csv']) );

end







% size of parfor loop
simulationSize = size (hIds ,1) %'1' here is for caculate the # of rows, which is the # of pairs


%% Setting arrays that storage the max growth rate. 'gv'->growth in vector format; 'gm'->growth in matrix format;
gv_hostMaxFirst                 = zeros(simulationSize ,1);
gv_guest_st_hostMaxFirst        = zeros(simulationSize ,1);


%% Define the directory path for INPUT: metabolic model 
cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL

dataUsed = '/CarveFungiEndo_rb';  % Define the location of input: metabolic model

% Define the input directory path
% C:\Users\lusa4312\OneDrive\ProkaryoteEndosymbiosis-main\CarveFungiEndo_rb;
filedirGEM = ['C:/Users/'  username  cloud_local  dataUsed ]                                                     



%% Finding parameters to use within loop
% Step1: Replaced EVAL with LOAD to load the host model, avoiding "transparency violation error"
host_model_path = fullfile(filedirGEM, ['model1.mat'])
ehmodel1 = load(host_model_path, 'metabolic_model');
ehmodel1 = ehmodel1.metabolic_model;


% Step2: Parameters used in function to create the GEM structure 
ne  = size(ehmodel1.S_ext      ,1) % ne  -> # of extracelular mapped metabolites
nic = size(ehmodel1.S_intC     ,1) % nic -> # of metabolites in citoplasm
nio = size(ehmodel1.S_intO     ,1) % nio -> # of metabolites in other compartments
nu  = size(ehmodel1.S_unmapped ,1) % nu  -> # of extracelular unmapped metabolites

% Step3: Cleaning
ehmodel1 =[];


%% Defining download directory
% Define where data is uploaded from: Onedrive (cloud) or local  
%cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL                        %% COMMENTED TO RUN ON THE CLUSTER

if strcmp (username, 'lusa4312')
    cloud_local = '/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
end

if strcmp (username, 'lucas')
    cloud_local = '/Desktop/OneDrive/ProkaryoteEndosymbiosis-main'; %CLOUD
end


% Define the output directory path
if strcmp (runningIn, 'desktop')
    %SavingPathName = 'C:\Users\lusa4312\OneDrive\ProkaryoteEndosymbiosis-main\growthResults_cfe\growthResultsTestD'; % alternative code for the one bellow 
    SavingPathName = ['C:\Users\' username cloud_local '\growthResults_cfe\growthResultsTestD']  %'D' stands for donation from 1st to 2nd
elseif strcmp (runningIn, 'hpc2n')   
    SavingPathName = '/pfs/proj/nobackup/fs/projnb10/hpc2n2023-112/lusa4312/Documents/ProkaryoteEndosymbiosis-main/growthResults_cfe/growthResultsTestD';
end




%% Initialize a parallel pool
parpool;                                                                                                                



parfor i = 1:simulationSize
    
    %% compute growth rates for host in own environs
    ind1 = hIds (i); % integer for a metabolic model for host
    ind2 = gIds(i); % integer for a metabolic model for endo
    Sync_Fitness = vec_P (i)
    
    % Load the host model
    % Change 1: Replaced EVAL with LOAD to load the host model, avoiding "transparency violation error"
    host_model_path = fullfile(filedirGEM, ['model', num2str(ind1), '.mat']);
    ehmodel1 = load(host_model_path, 'metabolic_model');
    ehmodel1 = ehmodel1.metabolic_model;
    

    % Load the guest (endo) model
    % Change 2: Replaced EVAL with LOAD to load the endo model, avoiding "transparency violation error"
    endo_model_path = fullfile(filedirGEM, ['model', num2str(ind2), '.mat']);
    ehmodel2 = load(endo_model_path, 'metabolic_model');
    ehmodel2 = ehmodel2.metabolic_model;


    %% Prefix used in saved name:
    hostMaxFirst_prefix = ['hID' num2str(ind1) '_gID'  num2str(ind2) '_wH1_of_H1G2_viab'  minGrowthStr  '_cfe_default_rhs_step_Sync']
    guestMaxSec_prefix  = ['hID' num2str(ind1) '_gID'  num2str(ind2) '_wG2_of_H1G2_viab'  minGrowthStr  '_cfe_default_rhs_step_Sync']



    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAXIMIXING THE HOST FIRST, THEN THE GUEST
    % Create endo model in two steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Step 1 of resultpair12: Maximize growth of host (M1), given constrain that 
    %                         [growth of guest (M2)] > 0  

    pairmodel = create_pair_step1_ofMHF_cfe_viab(ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu  ,minGrowth)
    gurobi_result_pair12_step1_ofMHF_firstM1 = gurobi(pairmodel,params);

    % check if feasible
    if ~strcmp(gurobi_result_pair12_step1_ofMHF_firstM1.status,'OPTIMAL')
        result_pair12_step1_ofMHF_firstM1 = 0;                           
    else
        result_pair12_step1_ofMHF_firstM1 = abs(gurobi_result_pair12_step1_ofMHF_firstM1.objval);
    end
    
        
    %% Step 2 of resultpair12: Maximize growth of guest (M2), using the maximized growth
    %                          of host (M1) as constrain

    gurobi_result_pair12_step2_ofMHF_firstM1 = struct(); % If result_pair12_step1_ofMHF_firstM1 is not greater than 0, 
                                                         % gurobi_result_pair12_step2_ofMHF_firstM1 is never set, 
                                                         % which would cause an error when you attempt to use 
                                                         % gurobi_result_pair12_step2_ofMHF_firstM1 later in the code. 
                                                         % One way to handle this is to initialize (gurobi_result_pair12_step2_ofMHF_firstM1 = struct()) 
                                                         % the variable with default values before your if conditions.

    if result_pair12_step1_ofMHF_firstM1 > thres_min_growth
        pairmodel = create_pair_step2_ofMHF_cfe (ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,result_pair12_step1_ofMHF_firstM1);
        gurobi_result_pair12_step2_ofMHF_firstM1 = gurobi(pairmodel,params);

        % check if feasible
        if ~strcmp(gurobi_result_pair12_step2_ofMHF_firstM1.status,'OPTIMAL')
            result_pair12_step2_ofMHF_firstM1 = 0;
        else
            result_pair12_step2_ofMHF_firstM1 = abs(gurobi_result_pair12_step2_ofMHF_firstM1.objval);
        end
    else
        result_pair12_step2_ofMHF_firstM1 = 0;
    end
       
      
    %% Calculate the number of intervals needed to reach or exceed X
    stepSize = 0.001  % increment value, this is how much the host invests to the guest's growth
    nIntervalsH = ceil(result_pair12_step1_ofMHF_firstM1 / stepSize);
    nIntervalsG = ceil(result_pair12_step2_ofMHF_firstM1 / stepSize);
    nIntervals  = max(nIntervalsH,nIntervalsG)

    % The estimated vector size is n intervals + 1 for the starting point at 0
    vecSize = nIntervals + 1;

    % preallocate space 
    gv_hostMaxFirst          = zeros(vecSize, 1);
    gv_guest_st_hostMaxFirst = zeros(vecSize, 1);

    % save wH & wG at t=1
    nRun = 1;
    maxRun = 30;
    gv_hostMaxFirst(nRun)          = result_pair12_step1_ofMHF_firstM1;
    gv_guest_st_hostMaxFirst(nRun) = result_pair12_step2_ofMHF_firstM1;


    % criteria to stop simulation
    stopThres = 0.001;

while (stopThres < result_pair12_step1_ofMHF_firstM1) && (stopThres < result_pair12_step2_ofMHF_firstM1) 

    nRun = nRun + 1

    % calculate at t>1
    %% Step 1 of resultpair12: Maximize growth of host (M1), given constrain that 
    %                        [growth of guest (M2)] > result_pair12_step2_ofMHF_firstM1 + stepSize  
    gurobi_result_pair12_step1_ofMHF_firstM1 = struct();
    result_pair12_step1_ofMHF_firstM1 = []; %cleaning to avoid coping past values by mistake.
    pairmodel = create_pair_step1_ofMHF_cfe_viab_step_Sync(ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,minGrowth ,result_pair12_step2_ofMHF_firstM1 ,stepSize  ,Sync_Fitness);
    gurobi_result_pair12_step1_ofMHF_firstM1 = gurobi(pairmodel,params);                                   

    % check if feasible
    if ~strcmp(gurobi_result_pair12_step1_ofMHF_firstM1.status,'OPTIMAL')
        result_pair12_step1_ofMHF_firstM1 = 0;                           
    else
        result_pair12_step1_ofMHF_firstM1 = abs(gurobi_result_pair12_step1_ofMHF_firstM1.objval);
    end

    %% Step 2 of resultpair12: Maximize growth of guest (M2), using the maximized growth
    %                          of host (M1) > result_pair12_step1_ofMHF_firstM1 + stepSize
    gurobi_result_pair12_step2_ofMHF_firstM1 = struct(); % If result_pair12_step1_ofMHF_firstM1 is not greater than 0, 
                                                         % gurobi_result_pair12_step2_ofMHF_firstM1 is never set, 
                                                         % which would cause an error when you attempt to use 
                                                         % gurobi_result_pair12_step2_ofMHF_firstM1 later in the code. 
                                                         % One way to handle this is to initialize (gurobi_result_pair12_step2_ofMHF_firstM1 = struct()) 
                                                         % the variable with default values before your if conditions.
    if result_pair12_step1_ofMHF_firstM1 > thres_min_growth
        result_pair12_step2_ofMHF_firstM1 = []; %cleaning to avoid coping past values by mistake.
        pairmodel = create_pair_step2_ofMHF_cfe (ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu ,result_pair12_step1_ofMHF_firstM1);   
                    
        gurobi_result_pair12_step2_ofMHF_firstM1 = gurobi(pairmodel,params);

        % check if feasible
        if ~strcmp(gurobi_result_pair12_step2_ofMHF_firstM1.status,'OPTIMAL')
            result_pair12_step2_ofMHF_firstM1 = 0;
        else
            result_pair12_step2_ofMHF_firstM1 = abs(gurobi_result_pair12_step2_ofMHF_firstM1.objval);
        end
    else
        result_pair12_step2_ofMHF_firstM1 = 0;
    end
    
    %% Check if resizing is needed
    if nRun > length(gv_hostMaxFirst)
        gv_hostMaxFirst          = [gv_hostMaxFirst;          zeros(vecSize, 1)]; % Double the size when needed
        gv_guest_st_hostMaxFirst = [gv_guest_st_hostMaxFirst; zeros(vecSize, 1)]; % Double the size when needed
    end
    
    % Store the total
    gv_hostMaxFirst(nRun)          = result_pair12_step1_ofMHF_firstM1;
    gv_guest_st_hostMaxFirst(nRun) = result_pair12_step2_ofMHF_firstM1;
end                  % end from 'while' loop  

% Trim the excess zeros if the estimate is too large
gv_hostMaxFirst          = gv_hostMaxFirst(1:nRun);
gv_guest_st_hostMaxFirst = gv_guest_st_hostMaxFirst(1:nRun);

%%

if gv_hostMaxFirst(end) == 0 %% ADDED in 2024.08.15
    % Remove the last row
    gv_hostMaxFirst          = gv_hostMaxFirst(1:end-1);
    gv_guest_st_hostMaxFirst = gv_guest_st_hostMaxFirst(1:end-1);
end


    %% Cleanning by overwriting with an empty array
    ehmodel2     = [];
    pairmodel    = [];
    % resultpairhe = [];
    % resulth2env  = [];
    result_pair12_step1_ofMHF_firstM1 = []; 
    result_pair12_step2_ofMHF_firstM1 = []; 
    %result_pair12_step22_ofMHF_firstM1_e_cannotUseUnusedHostFlux = [];
    gurobi_result_pair12_step1_ofMHF_firstM1 = struct();   
    gurobi_result_pair12_step2_ofMHF_firstM1 = struct();
    gurobi_result_pair12_step22_ofMHF_firstM1 = struct();
    %zero_flux_h = [];




    %% Cleanning by overwriting with an empty array
     metabolic_model = [];
     ehmodel1        = [];
     resulth         = [];
     hostmodel       = [];
     pairmodel       = [];


     %% Saving data 
 
     writematrix(gv_hostMaxFirst             ,fullfile(SavingPathName ,[hostMaxFirst_prefix      sufix_database '.csv']));
     writematrix(gv_guest_st_hostMaxFirst    ,fullfile(SavingPathName ,[guestMaxSec_prefix       sufix_database '.csv']));


end                                                                        % end of parfor loop




% Close the parallel pool
if strcmp (runningIn, 'desktop')
      delete(gcp('nocreate'));                                
elseif strcmp (runningIn, 'hpc2n')   
      %delete(gcp('nocreate'));                             %%COMMENTED TO RUN ON THE CLUSTER;
end


time = toc (t0);
disp ('Sim. Finished')





