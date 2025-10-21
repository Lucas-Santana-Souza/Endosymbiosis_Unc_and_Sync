% 2023.06.29 - Author: Lucas Santana Souza 
% Aim: to calculate the: 
%
%% [growth rate of the Host maximized first, s.t. (guest's growth) > 0]
%   HOST          GUEST
% [ S_ext ] | [     0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [ S_ext ] | [     0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [   0   ] | [S_unmapped ] [=] [0]                                           -> REGION b3  
% [ S_int ] | [ S_ext2int ] [=] [0]                                           -> REGION b4 
% [   0   ] | [   S_int   ] [=] [0]                                           -> REGION b5  
% [   0   ] | [     1     ] [>] [0]                                           -> REGION for constrain in the growth rate when host is maximized first, s.t. (guest's growth) > 0

%% [growth rate of the Guest maximized second, s.t. (host max first)]
%   HOST          GUEST
% [ S_ext ] | [     0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [ S_ext ] | [     0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [   0   ] | [S_unmapped ] [=] [0]                                           -> REGION b3  
% [ S_int ] | [ S_ext2int ] [=] [0]                                           -> REGION b4 
% [   0   ] | [   S_int   ] [=] [0]                                           -> REGION b5  
% [   1   ] | [     0     ] [>]  [result_pairij_step1_ofMHF_firstMi - ErrorTolerance]   -> REGION  for constrain in the growth rate when guest is maximized, s.t. host max first

% 'S_ext' -> Host's compartment only contain external metabolites that can be mapped
% 'S_int' -> Host's/Guest's compartment only contain internal metabolites 
% 'S_ext2int' -> Not sure if this guest's compartment contain external metabolites that can be mapped and unmapped or just mapped
% 'S_unmapped' ->  guest's compartment contain external metabolites that are unmapped 



tic;

clear;
clc;
format compact


%% Gurobi settings  
addpath 'C:\gurobi1001\win64\matlab'
params = struct();
params.OutputFlag = 0;
params.FeasibilityTol=1e-9;

%% Defining upload directory
% Define the username                                                      % change dependening on computer
%username = 'lsant';
%username =  'lucas';
username = 'lusa4312';

% Define where data is uploaded from: Onedrive (cloud) or local  
%cloud_local = '/OneDrive - Umeå universitet/ProkaryoteEndosymbiosis-main'; %CLOUD
cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL

% Define the data used
%dataUsed = '/ext_int_models_CarveMe';
dataUsed = '/ext_int_models_Agora';

% Define the input directory path
%filedir='C:\Users\lusa4312\Documents\ProkaryoteEndosymbiosis-main\renamingAGORA'; % alternative code for the one bellow 
filedir = ['C:/Users/'  username  cloud_local  dataUsed ];

%% Get a list of all files in the directory
files = dir(fullfile(filedir, '*.mat'));

% Count the number of ".mat" files
num_mat_files = length(files);

ArraySize = 5 %num_mat_files;                                            %CRITICAL PLACE TO RUN

%% Setting arrays that storage the max growth rate
gv_ancestral_alone_nonSharedEnv = zeros (ArraySize ,1);
gm_hostMaxFirst = zeros (ArraySize ,ArraySize);
gm_guest_st_hostMaxFirst = zeros (ArraySize ,ArraySize);

%% Initialize a parallel pool
parpool;

%%

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

        %
        % MAXIMIXING THE HOST FIRST, THEN THE GUEST

        % Create endo model in two steps
        % Step 1 of resultpair12: Maximize growth of M1, given constrain that 
        %                        (growth of M2) > 0  
       
        %clear pairmodel
        pairmodel = create_pair_step1_ofMHF_uei (ehmodel1 ,ehmodel2 ,ne ,ni, nu);
        result_pair12_step1_ofMHF_firstM1 = gurobi(pairmodel,params);

        % check if feasible
        if ~strcmp(result_pair12_step1_ofMHF_firstM1.status,'OPTIMAL')
            result_pair12_step1_ofMHF_firstM1 = 0;                           
        else
            result_pair12_step1_ofMHF_firstM1 = abs(result_pair12_step1_ofMHF_firstM1.objval);
        end
        
        
        % Step 2 of resultpair12: Maximize growth of M2, using the maximized growth
        %                         of M1 as constrain

        %clear pairmodel %cannot use "clear" inside parfor loop

        if result_pair12_step1_ofMHF_firstM1 > 0
            pairmodel = create_pair_step2_ofMHF_uei (ehmodel1 ,ehmodel2 ,ne ,ni, nu ,result_pair12_step1_ofMHF_firstM1);
            result_pair12_step2_ofMHF_firstM1 = gurobi(pairmodel,params);
    
            % check if feasible
            if ~strcmp(result_pair12_step2_ofMHF_firstM1.status,'OPTIMAL')
                result_pair12_step2_ofMHF_firstM1 = 0;
            else
                result_pair12_step2_ofMHF_firstM1 = abs(result_pair12_step2_ofMHF_firstM1.objval);
            end
        else
            result_pair12_step2_ofMHF_firstM1 = 0;
        end
        
        
        % stores the max growth rates
        gm_hostMaxFirst (i,j) = result_pair12_step1_ofMHF_firstM1;
        gm_guest_st_hostMaxFirst (i,j) = result_pair12_step2_ofMHF_firstM1;
        
        % Cleanning by overwriting with an empty array
        ehmodel2     = [];
        pairmodel    = [];
        % resultpairhe = [];
        % resulth2env  = [];
        result_pair12_step1_ofMHF_firstM1 = []; 
        result_pair12_step2_ofMHF_firstM1 = []; 

    end

    % Cleanning by overwriting with an empty array
     metabolic_model = [];
     ehmodel1        = [];
     resulth         = [];
     hostmodel       = [];
     pairmodel       = [];

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


%% Saving data based on whether the collection is 'CarveMe' or 'AGORA'
if strcmp (dataUsed, '/ext_int_models_CarveMe')
    writematrix(gv_ancestral_alone_nonSharedEnv ,fullfile(SavingPathName, 'max_growth_ancestor_alone_nonSharedEnv_uei_CarveMe.csv'));
    writematrix(gm_hostMaxFirst                 ,fullfile(SavingPathName, 'max_growth_H_st_growth_G_is_larger0_uei_CarveMe.csv'));                 
    writematrix(gm_guest_st_hostMaxFirst        ,fullfile(SavingPathName, 'max_growth_G_st_growth_H_is_maxFirst_uei_CarveMe.csv'));
elseif strcmp (dataUsed, '/ext_int_models_Agora')   
    writematrix(gv_ancestral_alone_nonSharedEnv ,fullfile(SavingPathName, 'max_growth_ancestor_alone_nonSharedEnv_uei_AGORA.csv'));
    writematrix(gm_hostMaxFirst                 ,fullfile(SavingPathName, 'max_growth_H_st_growth_G_is_larger0_uei_AGORA.csv'));                 
    writematrix(gm_guest_st_hostMaxFirst        ,fullfile(SavingPathName, 'max_growth_G_st_growth_H_is_maxFirst_uei_AGORA.csv'));
end

toc

