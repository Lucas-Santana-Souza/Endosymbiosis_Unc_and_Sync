% 2024.05.20 - Author: Lucas Santana Souza
% Aim: to create the following model
%% Structure of the holobiont matrix (created in the 'create_pair_step1_ofMHF_uei' function):
%    HOST          GUEST
% [ S_ext ] | [       0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )] -> REGION b1
% [ S_ext ] | [       0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )] -> REGION b2  
% [   0   ] | [Ng*S_unmapped] [=] [0]                                            -> REGION b3  
% [ S_int ] | [Ng*S_ext2int ] [=] [0]                                            -> REGION b4 
% [   0   ] | [   Ng*S_int  ] [=] [0]                                            -> REGION b5  
% [   0   ] | [       1     ] [>] [minGrowth]                                    -> REGION b6 for constrain in the growth rate when host is maximized first (h1Endo), s.t. (guest's growth) > minGrowth
% [   0   ] | [       1     ] [>] [result_pair12_step2_ofMHF_firstM1 + stepSize] -> REGION b7 for the constrain: (g2Endo's growth) > (result_pair12_step2_ofMHF_firstM1 + stepSize)
%
%                            ehmodel1 -> microbe that is the host
%                            ehmodel2 -> microbe that is the guest
%                             'S_ext' -> Host's compartment only contain external metabolites that can be mapped
%                             'S_int' -> Host's/Guest's compartment only contain internal metabolites 
%                         'S_ext2int' -> Not sure if this guest's compartment contain external metabolites that can be mapped and unmapped or just mapped
%                        'S_unmapped' -> guest's compartment contain external metabolites that are unmapped 
%                                  ne -> # of extracelular mapped metabolites
%                                  ni -> # of intracelullar metabolites
%                                  nu -> # of extracelular unmapped metabolites
%                                'Ng' -> # of guests
%                          'stepSize' -> increment for G's growth
% 'result_pair12_step2_ofMHF_firstM1' -> guest max 2nd in endosymbiosis

function endomodel= create_pair_step1_ofMHF_uei_Ng_viab_step(ehmodel1 ,ehmodel2 ,ne ,ni ,nu ,Ng ,minGrowth ,result_pair12_step2_ofMHF_firstM1 ,stepSize)


%% find biomass index
bmih=ehmodel1.bmi;
bmie=ehmodel2.bmi;

%% construct host S matrices for E and C compartments

%% Create .obj field, which determines what is maximized

nrh = size(ehmodel1.lb,1); %number of reactions a host has
nre = size(ehmodel2.lb,1); %number of reactions a guest has
f   = zeros(nrh+nre,1);
f(bmih) = -1; % host's growth rate is maximized
endomodel.obj = f;

%% Defining metabolite regions to be used at the combined stochiometric matrix (A) and combined rhs (.rhs)
b1 = [1:ne];             % row indices corresponding to A for extracellular lb
b2 = [1:ne] + ne;        % row indices corresponding to A for extracellular ub
b3 = [1:nu] + 2*ne;      % row indices corresponding to A for unmmapped extracellular ub
b4 = [1:ni] + 2*ne + nu;      % row indices corresponding to A for host's intracellular and guest's extracellular 
b5 = [1:ni] + ni + 2*ne + nu; % row indices corresponding to A for guest's intracellular 
b6 = 1 + (2*ne + nu + 2*ni); 
b7 = 1 + (2*ne + nu + 2*ni +1);

total_num_rows = 2*ne + nu + 2*ni + 2;

%% Pre-empty the combined stochiometric matrix (A)
endomat = sparse(total_num_rows ,nrh +nre);

%% Filling the combined stochiometric matrix (A)
%   HOST          GUEST
% [ S_ext ] | [     0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )] -> REGION b1
% [ S_ext ] | [     0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )] -> REGION b2  
% [   0   ] | [S_unmapped ] [=] [0]                                            -> REGION b3  
% [ S_int ] | [ S_ext2int ] [=] [0]                                            -> REGION b4 
% [   0   ] | [   S_int   ] [=] [0]                                            -> REGION b5  
% [   0   ] | [     1     ] [>] [minGrowth]                                    -> REGION b6 for constrain in the growth rate when host is maximized first (h1Endo), s.t. (guest's growth) > minGrowth
% [   0   ] | [     1     ] [>] [result_pair12_step2_ofMHF_firstM1 + stepSize] -> REGION b7 for the constrain: g2Endo's growth > result_pair12_step2_ofMHF_firstM1 + stepSize

% REGION b1: e compartment lower bound
endomat(b1 ,1:nrh) = ehmodel1.S_ext;
% endomat(b1 ,nrh+1:nre+nrh)= 0;        

% REGION b2: e compartment upper bound
endomat(b2 ,1:nrh) = ehmodel1.S_ext;
%endomat(b2 ,nrh+1:nre+nrh) = 0;        

% REGION b3:
%endomat(b3 ,1:nrh) = 0;                
endomat(b3 ,nrh+1:nre+nrh) = Ng*(ehmodel2.S_unmapped);

% REGION b4: host's intra. and guest's extrac. region 
endomat(b4 ,1:nrh) = ehmodel1.S_int;
endomat(b4 ,nrh+1:nre+nrh) = Ng*(ehmodel2.S_ext2int);

% REGION b5: c compartment endo
%endomat(b5 ,1:nrh) = 0;                
endomat(b5 ,nrh+1:nre+nrh) = Ng*(ehmodel2.S_int);




% REGION b6: for constrain in the growth rate 
endomat(b6 ,bmih)     = 0;                        
endomat(b6 ,bmie+nrh) = 1;                        

% REGION b7: for constrain g2Endo > result_pair12_step2_ofMHF_firstM1 + stepSize           
endomat(b7 ,bmih)     = 0;                                        
endomat(b7 ,bmie+nrh) = 1;   

% The combined stochiometric matrix (A) must be sparse to run in gurobi
endomodel.A = sparse(endomat);

%% Create .rhs field to use as condition for the combined stochiometric matrix (A) 

% Pre-empty field .rhs
endomodel.rhs = zeros(total_num_rows ,1);

% filing the field .rhs
endomodel.rhs(b1) = ehmodel1.rhs_ext_lb + ehmodel2.rhs_ext_lb;
endomodel.rhs(b2) = ehmodel1.rhs_ext_ub + ehmodel2.rhs_ext_ub;
%endomodel.rhs(b3) = zeros(nu, 1); %
%endomodel.rhs(b4) = zeros(ni, 1); % this makes the (host's intrac rhs) + (guest's extrac rhs) = 0
%endomodel.rhs(b5) = zeros(ni, 1); % this makes the (host's intrac rhs) + (guest's intrac rhs) = 0 
endomodel.rhs(b6) = minGrowth;             % modified here from 'endomodel.rhs(end)=0' (2024.05.12)--> 'This changes depending on who is max second or if the (Host's growth) = (Guest's growth)'
endomodel.rhs(b7) = result_pair12_step2_ofMHF_firstM1 + stepSize; %(Endo G's growth) > (result_pair12_step2_ofMHF_firstM1 + stepSize) addeed in (2024.05.20)


%% Create .sense field which defines the condition btw the matrix A and .rhs
endomodel.sense = [repmat('>' ,ne ,1); ... % -> REGION b1
                   repmat('<' ,ne ,1); ... % -> REGION b2
                   repmat('=' ,nu ,1); ... % -> REGION b3
                   repmat('=' ,ni ,1); ... % -> REGION b4
                   repmat('=' ,ni ,1); ... % -> REGION b5
                   repmat('>' ,1  ,1); ... % -> REGION b6 --> This is '=' if the (Host's growth) = (Guest's growth) --> modified in (2024.05.12)
                   repmat('>' ,1  ,1)];    % -> REGION b7 --> This is '>' as (Endo G's growth) = (result_pair12_step2_ofMHF_firstM1 + stepSize)   --> addeed in (2024.05.20)


%% Update fields for fluxes .lb & .ub
endomodel.lb=[ehmodel1.lb;ehmodel2.lb];                                                    
endomodel.ub=[ehmodel1.ub;ehmodel2.ub];

%% some models have an upper limit on biomass
endomodel.lb(bmih)=0;
endomodel.ub(bmih)=1000;
endomodel.lb(bmie+nrh)=0;
endomodel.ub(bmie+nrh)=1000;
end
