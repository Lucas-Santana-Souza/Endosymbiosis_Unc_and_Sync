% 2026.05.12 - Author: Lucas Santana Souza
%
% Aim: to create the following model
%% Structure of matrix for: Ancestral (here I will call HOST for illustration) Alone in Shared environment :
%   HOST         
% [ S_ext ]  [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [ S_ext ]  [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [ S_intC ] [=] [0]                                           -> REGION b3 
% [ S_intO ] [=] [0]                                           -> REGION b4 
%
%
%     ehmodel1 -> microbe that is the host
%     ehmodel2 -> microbe that is the guest
%      'S_ext' -> Host's compartment only contain external metabolites that can be mapped
%      'S_intC' -> Host's/Guest's only contain internal metabolites from the citoplasmatic compartment 
%      'S_intO' -> Host's/Guest's only contain internal metabolites from the other compartments (mitochodria, nucleos …)
%  'S_ext2int' -> Not sure if this guest's compartment contain external metabolites that can be mapped and unmapped or just mapped
% 'S_unmapped' ->  guest's compartment contain external metabolites that are unmapped 
%           ne -> # of extracelular mapped metabolites
%           nic -> # of intracelullar metabolites from S_intC
%           nio -> # of intracelullar metabolites from S_intO
%           nu -> # of extracelular unmapped metabolites


function endomodel= create_AloneShared_cfe (ehmodel1 ,ehmodel2 ,ne ,nic ,nio)



%% find biomass index
bmih=ehmodel1.bmi;
%bmie=ehmodel2.bmi;

%% construct host S matrices for E and C compartments

%% Create .obj field, which determines what is maximized
nrh = size(ehmodel1.lb,1); %number of reactions a host has
%nre = size(ehmodel2.lb,1); %number of reactions a guest has

f   = zeros(nrh,1);
f(bmih) = -1; % host's growth rate is maximized
endomodel.obj = f;

%% Defining metabolite regions to be used at the combined stochiometric matrix (A) and combined rhs (.rhs)
b1 =  1:ne;                % row indices corresponding to A for extracellular lb
b2 = (1:ne)  + ne;         % row indices corresponding to A for extracellular ub
b3 = (1:nic) + 2*ne;       % row indices corresponding to A for intracellular  (citoplasm)
b4 = (1:nio) + 2*ne + nic; % row indices corresponding to A for intracellular  (other compartments)

total_num_rows = 2*ne + nic + nio;

total_num_cols = nrh ;

%% Pre-empty the combined stochiometric matrix (A)
endomat = sparse(total_num_rows ,total_num_cols);

%% Filling the combined stochiometric matrix (A)
%   HOST          GUEST
% [ S_ext ]  [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [ S_ext ]  [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [ S_intC ] [=] [0]                                           -> REGION b3 
% [ S_intO ] [=] [0]                                           -> REGION b4 

% REGION b1: e compartment lower bound
endomat(b1 ,1:nrh)        =     ehmodel1.S_ext;
%endomat(b1 ,nrh+1:nre+nrh)=    ehmodel2.S_ext;

% REGION b2: e compartment upper bound
endomat(b2 ,1:nrh)         =     ehmodel1.S_ext;
%endomat(b2 ,nrh+1:nre+nrh) =    ehmodel2.S_ext; 

% REGION b3: Exosymbiosis host's intra. region 
endomat(b3 ,1:nrh)         = ehmodel1.S_intC;
% endomat(b3 ,nrh+1:nre+nrh) = ehmodel2.S_ext2int;

% REGION b4: Exosymbiosis host's intra. region 
endomat(b4 ,1:nrh)         = ehmodel1.S_intO;
% endomat(b4 ,nrh+1:nre+nrh) = 0;


% The combined stochiometric matrix (A) must be sparse to run in gurobi
endomodel.A = sparse(endomat);

%% Create .rhs field to use as condition for the combined stochiometric matrix (A) 

% Pre-empty field .rhs
endomodel.rhs = zeros(total_num_rows ,1);

% filing the field .rhs
endomodel.rhs(b1) = ehmodel1.rhs_ext_lb + ehmodel2.rhs_ext_lb;
endomodel.rhs(b2) = ehmodel1.rhs_ext_ub + ehmodel2.rhs_ext_ub;
%endomodel.rhs(b3) = zeros(nic, 1); 
%endomodel.rhs(b4) = zeros(nio, 1);  
%endomodel.rhs(end)=0;             

%% Create .sense field which defines the condition btw the matrix A and .rhs
endomodel.sense = [repmat('>' ,ne  ,1); ... % -> REGION b1
                   repmat('<' ,ne  ,1); ... % -> REGION b2
                   repmat('=' ,nic ,1); ... % -> REGION b3
                   repmat('=' ,nio ,1)];    % -> REGION b4

%% Update fields for fluxes .lb & .ub
endomodel.lb=[ehmodel1.lb]; %Obs: only hosts are present, if guests were also present this would change
endomodel.ub=[ehmodel1.ub]; %Obs: only hosts are present, if guests were also present this would change

%% some models have an upper limit on biomass
endomodel.lb(bmih)=0;
endomodel.ub(bmih)=1000;
end
