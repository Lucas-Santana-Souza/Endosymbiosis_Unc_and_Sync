% 2026.05.15 - Author: Lucas Santana Souza
%  AIM: for CarveFungiEndo
%   HOST          GUEST
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
%
function endomodel = create_pair_step1_ofMHF_cfe_viab_step_Sync(ehmodel1 ,ehmodel2 ,ne ,nic ,nio ,nu  ,minGrowth ,result_pair12_step2_ofMHF_firstM1 ,stepSize ,Sync_Fitness)

ErrorTolerance = 0.001;

%% find biomass index
bmih = ehmodel1.bmi;
bmie = ehmodel2.bmi;

%% construct host S matrices for E and C compartments

%% Create .obj field, which determines what is maximized

nrh = size(ehmodel1.lb,1); %number of reactions a host has
nre = size(ehmodel2.lb,1); %number of reactions a guest has

f   = zeros(nrh+nre,1);
f(bmih) = -1; % host's growth rate is maximized
endomodel.obj = f;

%% Defining metabolite regions to be used at the combined stochiometric matrix (A) and combined rhs (.rhs)
b1 = (1:ne);                           % row indices corresponding to A for extracellular lb
b2 = (1:ne) + ne;                      % row indices corresponding to A for extracellular ub
b3 = (1:nu) + 2*ne;                    % row indices corresponding to A for unmmapped extracellular 
b4 = (1:nic) + 2*ne + nu;              % row indices corresponding to A for host's intracellular and guest's extracellular 
b5 = (1:nic) + nic + 2*ne + nu;        % row indices corresponding to A for guest's intracellular C
b6 = (1:nio) + 2*nic + 2*ne + nu;      % row indices corresponding to A for guest's intracellular O  
b7 = (1:nio) + nio +2*nic + 2*ne + nu; % row indices corresponding to A for host's  intracellular O
b8 =  1 + (2*nio +2*nic + 2*ne + nu);      % row indices corresponding to A for guest's growth > minGrowth
b9 =  1 + (1 + 2*nio +2*nic + 2*ne + nu);  % row indices corresponding to A for guest's growth constraint
b10 = 1 + (2 + 2*nio +2*nic + 2*ne + nu);  % row indices corresponding to A for host's growth constraint

total_num_rows = 2*ne + nu + 2*nic + 2*nio + 3;

%% Pre-empty the combined stochiometric matrix (A)
endomat = sparse(total_num_rows ,nrh +nre);

%% Filling the combined stochiometric matrix (A)
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

% REGION b1: e compartment lower bound
endomat(b1 ,1:nrh) = ehmodel1.S_ext;
% endomat(b1 ,nrh+1:nre+nrh)= 0;       

% REGION b2: e compartment upper bound
endomat(b2 ,1:nrh) = ehmodel1.S_ext;
%endomat(b2 ,nrh+1:nre+nrh) = 0;        

% REGION b3: guest's S_unmapped
%endomat(b3 ,1:nrh) = 0;                
endomat(b3 ,nrh+1:nre+nrh) = ehmodel2.S_unmapped;

% REGION b4: host's intra. and guest's extrac. region 
endomat(b4 ,1:nrh) = ehmodel1.S_intC;
endomat(b4 ,nrh+1:nre+nrh) = ehmodel2.S_ext2intC;

% REGION b5: guest's S_intC compartment 
%endomat(b5 ,1:nrh) = 0;                
endomat(b5 ,nrh+1:nre+nrh) = ehmodel2.S_intC;

% REGION b6: guest's S_intO compartment
%endomat(b6 ,1:nrh) = 0;                
endomat(b6 ,nrh+1:nre+nrh) = ehmodel2.S_intO;

% REGION b7: host's S_intO compartment
endomat(b7 ,1:nrh) = ehmodel1.S_intO;
%endomat(b7 ,nrh+1:nre+nrh) = 0;        

% REGION8 for  (guest's growth) > minGrowth
%endomat(b8 ,bmih)     =  0;
endomat(b8 ,bmie+nrh) = 1; % ensures they grow at same rate

% REGION b9: for constrain g2Endo > result_pair12_step2_ofMHF_firstM1 + stepSize           
%endomat(b9 ,bmih)     = 0;                                        
endomat(b9 ,bmie+nrh) = 1;   

% REGION b10: for constrain h1Endo > Sync_Fitness - ErrorTolerance           
endomat(b10 ,bmih)     = 1;                                        
%endomat(b10 ,bmie+nrh) = 0;   

% The combined stochiometric matrix (A) must be sparse to run in gurobi
endomodel.A = sparse(endomat);


%% Create .rhs field to use as condition for the combined stochiometric matrix (A) 
% Pre-empty field .rhs
endomodel.rhs = zeros(total_num_rows ,1);

% filing the field .rhs
endomodel.rhs(b1)  = ehmodel1.rhs_ext_lb + ehmodel2.rhs_ext_lb;
endomodel.rhs(b2)  = ehmodel1.rhs_ext_ub + ehmodel2.rhs_ext_ub;
%endomodel.rhs(b3) = zeros(nu, 1); %
%endomodel.rhs(b4) = zeros(nic, 1); %
%endomodel.rhs(b5) = zeros(nic, 1); %
%endomodel.rhs(b6) = zeros(nio, 1); %
%endomodel.rhs(b7) = zeros(nio, 1); %
endomodel.rhs(b8)  = minGrowth; %
endomodel.rhs(b9)  = result_pair12_step2_ofMHF_firstM1 + stepSize; %(Endo G's growth) > (result_pair12_step2_ofMHF_firstM1 + stepSize)  
endomodel.rhs(b10) = Sync_Fitness - ErrorTolerance;                %(Endo H's growth) > (Sync_Fitness - ErrorTolerance)  


%% Create .sense field which defines the condition btw the matrix A and .rhs
endomodel.sense = [repmat('>' ,ne  ,1); ...  % -> REGION b1
                   repmat('<' ,ne  ,1); ...  % -> REGION b2
                   repmat('=' ,nu  ,1); ...  % -> REGION b3
                   repmat('=' ,nic ,1); ...  % -> REGION b4
                   repmat('=' ,nic ,1); ...  % -> REGION b5
                   repmat('=' ,nio ,1); ...  % -> REGION b6
                   repmat('=' ,nio ,1); ...  % -> REGION b7
                   repmat('>' ,1   ,1); ...  % -> REGION b8  This is '>' for the (Guest's growth) > minGrowth
                   repmat('>' ,1   ,1); ...  % -> REGION b9  This is '>' as (Endo G's growth) = (result_pair12_step2_ofMHF_firstM1 + stepSize)    
                   repmat('>' ,1   ,1)];     % -> REGION b10 This is '>' as (Endo H's growth) = (Sync_Fitness - ErrorTolerance)    

%% Update fields for fluxes .lb & .ub
endomodel.lb=[ehmodel1.lb;ehmodel2.lb];
endomodel.ub=[ehmodel1.ub;ehmodel2.ub];

% some models have an upper limit on biomass
endomodel.lb(bmih)=0;
endomodel.ub(bmih)=1000;
endomodel.lb(bmie+nrh)=0;
endomodel.ub(bmie+nrh)=1000;
end
