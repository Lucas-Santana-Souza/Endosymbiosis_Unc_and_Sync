function [time ,num_mat_files ,gv_ancestral_alone_nonSharedEnv,gv_hostMaxFirst ,gv_guest_st_hostMaxFirst] = f_MaxEndoHostFirst_then_guest_uei_Ng_withoutAloneS_viab_step(Ngs_input)
% 2023.05.20 - Author: Lucas Santana Souza 
% Aim: to calculate growth rates for hosts having increasing investment (stepSize) in guest's growth
%      which will be used to calculate trande-off curves.   
% 
%% create_pair_step1_ofMHF_uei_Ng_viab(ehmodel1 ,ehmodel2 ,ne ,ni ,nu ,Ng ,minGrowth); 
% [growth rate of the Host maximized first, s.t. (guest's growth) > minGrowth]
%   HOST          GUEST
% [ S_ext ] | [       0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [ S_ext ] | [       0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [   0   ] | [Ng*S_unmapped] [=] [0]                                           -> REGION b3  
% [ S_int ] | [Ng*S_ext2int ] [=] [0]                                           -> REGION b4 
% [   0   ] | [   Ng*S_int  ] [=] [0]                                           -> REGION b5 
% [   0   ] | [      1      ] [>] [minGrowth]                                   -> REGION b6 for constrain in the growth rate when host is maximized first (h1Endo), s.t. (guest's growth) > minGrowth
%
%% create_pair_step2_ofMHF_uei_Ng(ehmodel1 ,ehmodel2 ,ne ,ni, nu ,result_pair12_step1_ofMHF_firstM1 ,Ng); 
% [growth rate of the Guest maximized second, s.t. (host max first)]
%   HOST          GUEST
% [ S_ext ] | [       0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]        -> REGION b1
% [ S_ext ] | [       0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]        -> REGION b2  
% [   0   ] | [Ng*S_unmapped] [=] [0]                                                   -> REGION b3  
% [ S_int ] | [Ng*S_ext2int ] [=] [0]                                                   -> REGION b4 
% [   0   ] | [   Ng*S_int  ] [=] [0]                                                   -> REGION b5  
% [   1   ] | [       0     ] [>] [result_pairij_step1_ofMHF_firstMi - ErrorTolerance]  -> REGION b6 for constrain: guest's maximized growth rate is s.t. host's max first
%
%% create_pair_step1_ofMHF_uei_Ng_viab_step(ehmodel1 ,ehmodel2 ,ne ,ni ,nu ,Ng ,minGrowth ,result_pair12_step2_ofMHF_firstM1 ,stepSize); 
% [growth rate of the Host maximized first, s.t. (guest's growth) > minGrowth]
%   HOST          GUEST
% [ S_ext ] | [       0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )] -> REGION b1
% [ S_ext ] | [       0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )] -> REGION b2  
% [   0   ] | [Ng*S_unmapped] [=] [0]                                            -> REGION b3  
% [ S_int ] | [Ng*S_ext2int ] [=] [0]                                            -> REGION b4 
% [   0   ] | [   Ng*S_int  ] [=] [0]                                            -> REGION b5 
% [   0   ] | [      1      ] [>] [minGrowth]                                    -> REGION b6 for constrain in the growth rate when host is maximized first (h1Endo), s.t. (guest's growth) > minGrowth
% [   0   ] | [      1      ] [>] [result_pair12_step2_ofMHF_firstM1 + stepSize] -> REGION b7 for the constrain: g2Endo's growth > (result_pair12_step2_ofMHF_firstM1 + stepSize)
%
%
%
% endomodel.lb(zero_flux_h) = 0;
% endomodel.ub(zero_flux_h) = 0;
%
% 'S_ext' -> Host's compartment only contain external metabolites that can be mapped
% 'S_int' -> Host's/Guest's compartment only contain internal metabolites 
% 'S_ext2int' -> Not sure if this guest's compartment contain external metabolites that can be mapped and unmapped or just mapped
% 'S_unmapped' ->  guest's compartment contain external metabolites that are unmapped 
% 'zero_flux_h' -> indices of host's reactions that the host did not use when optimized first 
% 'Ng' -> # of guests
%


format compact
disp ('Sim. Started')
t0 = tic;

%% Defining where I am running
runningIn = 'desktop';
%runningIn = 'hpc2n';

%% Defyning if I am running a test simulation (a subset of metabolic models)
isThisATestSimulation = 'yes';
%isThisATestSimulation = 'no';

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

%% To be used in prefix of saved name:
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
%username = 'lsant';
%username =  'lucas';
username = 'lusa4312';

% Define the collection used
%collection = 'CarveMe'
collection = 'Agora'

%% DEFINE Database for upload growth rates already calculated
% Define where data is uploaded from: Onedrive (cloud) or local  
%cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL
cloud_local = '/OneDrive - Umeå universitet/ProkaryoteEndosymbiosis-main'; %CLOUD

if strcmp (collection, 'CarveMe')
   database = '/growthResults_uei/growthResultsCarveMe/';  % Define the location of input: growth of specified collection
   sufix_database = '_CarveMe';                            % Define the sufix of input file
elseif strcmp (collection, 'Agora')   
   database = '/growthResults_uei/growthResultsAgora/'; 
   sufix_database = '_AGORA';
end

%% Define the directory path for INPUT: growth rate in exosymbiosis
%dir_path = 'C:/Users/lusa4312/Documents/ProkaryoteEndosymbiosis-main/growthResults/growthResultsAgora/test/';
if strcmp (runningIn, 'desktop')
   dir_path  = ['C:/Users/'  username  cloud_local  database 'data_in_format_to_analyse']; 
elseif strcmp (runningIn, 'hpc2n')   
   dir_path = ['/pfs/proj/nobackup/fs/projnb10/hpc2n2023-112/lusa4312/Documents/ProkaryoteEndosymbiosis-main' database 'data_in_format_to_analyse'];
end












%% Upload hostID and guestID for cases in which: (h_Endo > h_Exo) && (g_Endo > g_Exo)

% Count the number of ".mat" files
if strcmp (isThisATestSimulation, 'yes')
      pairList =  [21 35; 566 6];  %[21 35];%                                         % Example values for AGORA when H1Endo>H1Exo & G2Endo>G2Exo
elseif strcmp (isThisATestSimulation, 'no')   
      input_file_name = ['hId_gId_H1G2_endo_viab' minGrowthStr '_Mutual.csv'];            
      hId_gId = fullfile(dir_path, input_file_name);   
      pairList = readmatrix(hId_gId);                                      % Data structure: 1st-column is hostID, 2nd-column is guestID  
end

% To avoid communication overhead (when using the index of the parfor loop): 
% by separating the data of IDs 
hIds = pairList (:,1);
gIds = pairList (:,2);

% size of parfor loop
num_mat_files = size (pairList,1) %'1' here is for caculate the # of rows, which is the # of pairs

%% Define the directory path for INPUT: metabolic model 
% cloud_local = '/OneDrive - Umeå universitet/ProkaryoteEndosymbiosis-main'; %CLOUD
cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL

if strcmp (collection, 'CarveMe')
   dataUsed = '/ext_int_models_CarveMe';  % Define the location of input: metabolic model
elseif strcmp (collection, 'Agora')   
   dataUsed = '/ext_int_models_Agora'; 
end

% Define the input directory path
if strcmp (runningIn, 'desktop')
      %filedir='C:\Users\lusa4312\Documents\ProkaryoteEndosymbiosis-main\renamingAGORA'; % alternative code for the one bellow 
      filedir = ['C:/Users/'  username  cloud_local  dataUsed ];                                                     
elseif strcmp (runningIn, 'hpc2n')   
      filedir = ['/pfs/proj/nobackup/fs/projnb10/hpc2n2023-112/'  username  cloud_local  dataUsed] ;%hpc2n;
end














%% Setting arrays that storage the max growth rate. 'gv'->growth in vector format; 'gm'->growth in matrix format;
gv_ancestral_alone_nonSharedEnv = zeros(1 ,1); 
gv_hostMaxFirst                 = zeros(1 ,1);
gv_guest_st_hostMaxFirst        = zeros(1 ,1);


%% Parameter initialization

Ngs = Ngs_input; %[10, 100, 1000]%2%1 % [1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10];  %Ng = 1; % Ng = Number of guests

[ne, ni, nu] = calc_ne_ni_nu (filedir);
% ne -> # of extracelular mapped metabolites
% ni -> # of intracelullar metabolites 
% nu -> # of extracelular unmapped metabolites

%% Initialize a parallel pool
if strcmp (runningIn, 'desktop')
      parpool;                                                                                                                
elseif strcmp (runningIn, 'hpc2n')   
      %parpool;                                                            %%COMMENTED TO RUN ON THE CLUSTER
end

for NG_id = 1:length(Ngs)

   Ng = Ngs (NG_id) 

parfor i = 1:num_mat_files
    
    %% compute growth rates for host in own environs
    ind1 = hIds (i); % integer for a metabolic model for host
    
    % Load the host model
    % Change 1: Replaced EVAL with LOAD to load the host model, avoiding "transparency violation error"
    host_model_path = fullfile(filedir, ['model', num2str(ind1), '.mat']);
    ehmodel1 = load(host_model_path, 'metabolic_model');
    ehmodel1 = ehmodel1.metabolic_model;
    








        
        ind2 = gIds(i); % integer for a metabolic model for endo
        
        % Load the guest (endo) model
        % Change 2: Replaced EVAL with LOAD to load the endo model, avoiding "transparency violation error"
        endo_model_path = fullfile(filedir, ['model', num2str(ind2), '.mat']);
        ehmodel2 = load(endo_model_path, 'metabolic_model');
        ehmodel2 = ehmodel2.metabolic_model;



        

        %% Prefix used in saved name:
        hostMaxFirst_prefix = ['hID' num2str(ind1) '_gID'  num2str(ind2) '_wH1_of_H1G2_viab'  minGrowthStr  '_uei_default_rhs_step_Ng']
        guestMaxSec_prefix  = ['hID' num2str(ind1) '_gID'  num2str(ind2) '_wG2_of_H1G2_viab'  minGrowthStr  '_uei_default_rhs_step_Ng']

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAXIMIXING THE HOST FIRST, THEN THE GUEST
        % Create endo model in two steps
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Step 1 of resultpair12: Maximize growth of host (M1), given constrain that 
        %                         [growth of guest (M2)] > 0  

        pairmodel = create_pair_step1_ofMHF_uei_Ng_viab(ehmodel1 ,ehmodel2 ,ne ,ni ,nu ,Ng ,minGrowth);  %******TO CHANGE HERE ********
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
            pairmodel = create_pair_step2_ofMHF_uei_Ng(ehmodel1 ,ehmodel2 ,ne ,ni, nu ,result_pair12_step1_ofMHF_firstM1 ,Ng);  %******TO CHANGE HERE ********
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

while (stopThres < result_pair12_step1_ofMHF_firstM1) && (stopThres < result_pair12_step2_ofMHF_firstM1) %&& (nRun < maxRun) % uncoment if i want to impose time constraint

    nRun = nRun + 1

    % calculate at t>1
        %% Step 1 of resultpair12: Maximize growth of host (M1), given constrain that 
        %                        [growth of guest (M2)] > result_pair12_step2_ofMHF_firstM1 + stepSize  
        gurobi_result_pair12_step1_ofMHF_firstM1 = struct();
        result_pair12_step1_ofMHF_firstM1 = []; %cleaning to avoid coping past values by mistake.
        pairmodel = create_pair_step1_ofMHF_uei_Ng_viab_step(ehmodel1 ,ehmodel2 ,ne ,ni ,nu ,Ng ,minGrowth ,result_pair12_step2_ofMHF_firstM1 ,stepSize);  %******TO CHANGE HERE ******** Error using f_MaxEndoHostFirst_then_guest_uei_Ng_withoutAloneS_viab_step Too many input arguments
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
            pairmodel = create_pair_step2_ofMHF_uei_Ng(ehmodel1 ,ehmodel2 ,ne ,ni, nu ,result_pair12_step1_ofMHF_firstM1 ,Ng);  %******TO CHANGE HERE ********
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
end                  % end from 'while (stopThres < result_pair12_step1_ofMHF_firstM1) && (stopThres < result_pair12_step2_ofMHF_firstM1) && (nRun < maxRun)'

% Trim the excess zeros if the estimate is too large
gv_hostMaxFirst          = gv_hostMaxFirst(1:nRun);
gv_guest_st_hostMaxFirst = gv_guest_st_hostMaxFirst(1:nRun);

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

 %   end                                                                   %TO DELETE for j loop


    %% Cleanning by overwriting with an empty array
     metabolic_model = [];
     ehmodel1        = [];
     resulth         = [];
     hostmodel       = [];
     pairmodel       = [];


     %% Defining download directory
     % Define where data is uploaded from: Onedrive (cloud) or local  
     %cloud_local = '/Documents/ProkaryoteEndosymbiosis-main'; %LOCAL                        %% COMMENTED TO RUN ON THE CLUSTER
     cloud_local = '/OneDrive - Umeå universitet/ProkaryoteEndosymbiosis-main'; %CLOUD       %% COMMENTED TO RUN ON THE CLUSTER
     
     % Define the output directory path
     if strcmp (runningIn, 'desktop')
      %SavingPathName = 'C:\Users\lusa4312\Documents\ProkaryoteEndosymbiosis-main\growthResults\growthResultsTest'; % alternative code for the one bellow 
      SavingPathName = ['C:\Users\' username cloud_local '\growthResults_uei\growthResultsTestD']; %'D' stands for donation from 1st to 2nd
     elseif strcmp (runningIn, 'hpc2n')   
      SavingPathName = '/pfs/proj/nobackup/fs/projnb10/hpc2n2023-112/lusa4312/Documents/ProkaryoteEndosymbiosis-main/growthResults_uei/growthResultsTestD';
     end

     %% Saving data based on whether the collection is 'CarveMe' or 'AGORA'
     if strcmp (dataUsed, '/ext_int_models_CarveMe')
         sufix_database = '_CarveMe';  
         %writematrix(gv_ancestral_alone_nonSharedEnv     ,fullfile(SavingPathName ,[aloneNonShared_prefix  num2str(Ng)  sufix_database '.csv']));
         writematrix(gv_hostMaxFirst                     ,fullfile(SavingPathName ,[hostMaxFirst_prefix    num2str(Ng)  sufix_database '.csv']));
         writematrix(gv_guest_st_hostMaxFirst            ,fullfile(SavingPathName ,[guestMaxSec_prefix     num2str(Ng)  sufix_database '.csv']));
     elseif strcmp (dataUsed, '/ext_int_models_Agora')   
         sufix_database = '_AGORA';
         %writematrix(gv_ancestral_alone_nonSharedEnv     ,fullfile(SavingPathName ,[aloneNonShared_prefix  num2str(Ng)  sufix_database '.csv']));
         writematrix(gv_hostMaxFirst                     ,fullfile(SavingPathName ,[hostMaxFirst_prefix    num2str(Ng)  sufix_database '.csv']));
         writematrix(gv_guest_st_hostMaxFirst            ,fullfile(SavingPathName ,[guestMaxSec_prefix     num2str(Ng)  sufix_database '.csv']));
     end


end                                                                        % end of parfor loop


end   % end from 'for NG_id = 1:length(Ngs)' line 188

% Close the parallel pool
if strcmp (runningIn, 'desktop')
      delete(gcp('nocreate'));                                
elseif strcmp (runningIn, 'hpc2n')   
      %delete(gcp('nocreate'));                             %%COMMENTED TO RUN ON THE CLUSTER;
end


time = toc (t0);
disp ('Sim. Finished')

%end



















