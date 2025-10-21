% 2023.06.22 - Author: Lucas Santana Souza 
% 
% Aim: save data for max growth in different conditions: 
%      [ancestor in sole env]; 
%      [pair with equal growth], 
%      [ancestors in shared env]
%
% The reason for changing the code:
% To separate external and internal compounds in different
% matrices, this way making easier to manipulate and shortening the matrix
% size.
%
%% Structure of the holobiont matrix (created in the 'create_holobiont_default_rhs.m' function):
%   HOST          GUEST
% [ S_ext ] | [     0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [ S_ext ] | [     0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [   0   ] | [S_unmapped ] [=] [0]                                           -> REGION b3  
% [ S_int ] | [ S_ext2int ] [=] [0]                                           -> REGION b4 
% [   0   ] | [   S_int   ] [=] [0]                                           -> REGION b5  
% [   1   ] | [    -1     ] [=] [0]   -> REGION for s.t (Host's growth) = (Guest's growth)
%
% 'S_ext' -> Host's compartment only contain external metabolites that can be mapped
% 'S_int' -> Host's compartment only contain internal metabolites 
% 'S_ext2int' -> Not sure if this guest's compartment contain external metabolites that can be mapped and unmapped or just mapped
% 'S_unmapped' ->  guest's compartment contain external metabolites that are unmapped 



tic;

clear
clc
format compact

%% Gurobi settings  
addpath 'C:\gurobi1001\win64\matlab'
params = struct();
params.OutputFlag = 0;
params.FeasibilityTol= 1e-9;
params.OptimalityTol = 1e-9; % Default value in gurobi is: 1e-6

%% Defining upload directory
% Define the username                                                      % change dependening on computer
%username = 'lsant';
%username =  'lucas';
username = 'lusa4312';

% Define where data is uploaded from: Onedrive (cloud) or local  
%cloud_local = '/OneDrive - Umeå universitet/ProkaryoteEndosymbiosis-main'; %CLOUD
cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL

% Define the data used
dataUsed = '/ext_int_models_CarveMe';
%dataUsed = '/ext_int_models_Agora';

% Define the input directory path
%filedir='C:\Users\lusa4312\Documents\ProkaryoteEndosymbiosis-main\renamingAGORA'; % alternative code for the one bellow 
filedir = ['C:/Users/'  username  cloud_local  dataUsed ];


%% Get a list of all files in the directory
files = dir(fullfile(filedir, '*.mat'));

% Count the number of ".mat" files
num_mat_files = length(files);

ArraySize = 5 %num_mat_files;  %  10;%                                         %CRITICAL PLACE TO RUN

%% Setting arrays that storage the max growth rate
gv_ancestral_alone_nonSharedEnv = zeros (ArraySize ,1);
gm_pair_st_equal_growth = zeros (ArraySize ,ArraySize);
gm_ancestral_alone_SharedEnv = zeros (ArraySize ,ArraySize);


%% Initialize a parallel pool
parpool;

%% Calculation
parfor i = 1:ArraySize
    
    
    % compute growth rates for host in own environs
    ind1 = i; % integer for a metabolic model for host
    
    % Load the host model
    % Change 1: Replaced EVAL with LOAD to load the host model, avoiding "transparency violation error"
    host_model_path = fullfile(filedir, ['model', num2str(ind1), '.mat']);
    ehmodel1 = load(host_model_path, 'metabolic_model');
    ehmodel1 = ehmodel1.metabolic_model;
    
    % Outputs host's growth rate s.t. growing alone in its own env (non-shared env)
    [resulth,hostmodel] = run_gh_model_ext_int(ehmodel1); 
    
    % Stores the max growth rates of ancestors growing alone in nonShared env
    gv_ancestral_alone_nonSharedEnv(i,1)= abs(resulth.objval); 
    
    % Parameters used in the function 'create_holobiont_default_rhs'
    ne = size(ehmodel1.S_ext      ,1); % ne -> # of extracelular mapped metabolites
    ni = size(ehmodel1.S_int      ,1); % ni -> # of intracelullar metabolites 
    nu = size(ehmodel1.S_unmapped ,1); % nu -> # of extracelular unmapped metabolites

    for j = 1:ArraySize
        
        ind2 = j; % integer for a metabolic model for endo
        
        % Load the guest (endo) model
        % Change 2: Replaced EVAL with LOAD to load the endo model, avoiding "transparency violation error"
        endo_model_path = fullfile(filedir, ['model', num2str(ind2), '.mat']);
        ehmodel2 = load(endo_model_path, 'metabolic_model');
        ehmodel2 = ehmodel2.metabolic_model;

        % Create holobiont model s.t. (host's growth rate) = (guest's growth rate) 
        pairmodel = create_holobiont_default_rhs (ehmodel1,ehmodel2 ,ne ,ni, nu);
        
        % Outputs host's growth rate s.t. (host's growth rate) = (guest's growth rate)
        resultpairhe = gurobi(pairmodel,params);
        % Check if feasible
        if ~strcmp(resultpairhe.status,'OPTIMAL')
            resultpairhe=0;
        else
            resultpairhe=abs(resultpairhe.objval);
        end

        % Grow host alone in Shared Environment
        hostmodel.rhs_ext_lb = pairmodel.rhs(1 : ne);      % make (host's rhs_ext_lb) = (holo's rhs corresponding ext lb)
        hostmodel.rhs_ext_ub = pairmodel.rhs(ne+1 : 2*ne); % make (host's rhs_ext_ub) = (holo's rhs corresponding ext ub)
        
        % Outputs host's growth rate s.t. growing alone in shared environment
        resulth2env = run_gh_model_ext_int(hostmodel); 

        % Stores the max growth rates
        gm_pair_st_equal_growth (i,j) = resultpairhe;
        gm_ancestral_alone_SharedEnv (i,j) = abs(resulth2env.objval);

        % Clean structs
        ehmodel2     = [];
        pairmodel    = [];
        resultpairhe = [];
        resulth2env  = [];

    end

    % Clean structs
    metabolic_model = [];
    ehmodel1        = [];
    resulth         = [];
    hostmodel       = [];

end

% Close the parallel pool
delete(gcp('nocreate'));

%% Defining download directory

% Define where data is uploaded from: Onedrive (cloud) or local  
%cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL
cloud_local = '/OneDrive - Umeå universitet/ProkaryoteEndosymbiosis-main'; %CLOUD

% Define the output directory path
%SavingPathName = 'C:\Users\lusa4312\Documents\ProkaryoteEndosymbiosis-main\growthResults\growthResultsTest'; % alternative code for the one bellow 
SavingPathName = ['C:\Users\' username cloud_local '\growthResults\growthResultsTest'];

%% Saving data

if strcmp (dataUsed, '/ext_int_models_CarveMe')
    writematrix(gv_ancestral_alone_nonSharedEnv, fullfile( SavingPathName, 'max_growth_ancestor_alone_nonSharedEnv_default_rhs_parallel_uei_CarveMe.csv'));
    writematrix(gm_pair_st_equal_growth,         fullfile( SavingPathName, 'max_growth_pair_st_equal_growth_default_rhs_parallel_uei_CarveMe.csv'));
    writematrix(gm_ancestral_alone_SharedEnv,    fullfile( SavingPathName, 'max_growth_ancestor_alone_SharedEnv_default_rhs_parallel_uei_CarveMe.csv'));
elseif strcmp (dataUsed, '/ext_int_models_Agora')   
    writematrix(gv_ancestral_alone_nonSharedEnv, fullfile( SavingPathName, 'max_growth_ancestor_alone_nonSharedEnv_default_rhs_parallel_uei_Agora.csv'));
    writematrix(gm_pair_st_equal_growth,         fullfile( SavingPathName, 'max_growth_pair_st_equal_growth_default_rhs_parallel_uei_Agora.csv'));
    writematrix(gm_ancestral_alone_SharedEnv,    fullfile( SavingPathName, 'max_growth_ancestor_alone_SharedEnv_default_rhs_parallel_uei_Agora.csv'));
end



toc
