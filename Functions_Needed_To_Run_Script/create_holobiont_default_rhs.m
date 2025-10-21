% 2023.06.11 - Author: Lucas Santana Souza
% Aim: to create the following model
%   HOST          GUEST
% [ S_ext ] | [     0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [ S_ext ] | [     0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [   0   ] | [S_unmapped ] [=] [0]                                           -> REGION b3  
% [ S_int ] | [ S_ext2int ] [=] [0]                                           -> REGION b4 
% [   0   ] | [   S_int   ] [=] [0]                                           -> REGION b5  
% [   1   ] | [    -1     ] [=] [0]   -> REGION for s.t (Host's growth) = (Guest's growth)
%
% hemodel1 -> microbe that is the host
% hemodel2 -> microbe that is the guest
%       ne -> # of extracelular mapped metabolites
%       ni -> # of intracelullar metabolites
%       nu -> # of extracelular unmapped metabolites



function endomodel = create_holobiont_default_rhs (hemodel1 ,hemodel2 ,ne ,ni ,nu)

  
% find biomass index
bmih = hemodel1.bmi;
bmie = hemodel2.bmi;


%% construct host S matrices for E and C compartments


%% Create .obj field, which determines what is maximized

nrh = size(hemodel1.lb,1); %number of reactions a host has
nre = size(hemodel2.lb,1); %number of reactions a guest has
f   = zeros(nrh+nre,1);
f(bmih) = -1; % host's growth rate is maximized
endomodel.obj = f;

%% Defining metabolite regions to be used at the combined stochiometric matrix (A) and combined rhs (.rhs)
b1 = [1:ne];             % row indices corresponding to S for extracellular lb
b2 = [1:ne] + ne;        % row indices corresponding to S for extracellular ub
b3 = [1:nu] + 2*ne;
b4 = [1:ni] + 2*ne + nu;      % row indices corresponding to S for host's intracellular and guest's extracellular 
b5 = [1:ni] + ni + 2*ne + nu; % row indices corresponding to S for guest's intracellular 

total_num_rows = 2*ne + nu + 2*ni + 1;

%% Pre-empty the combined stochiometric matrix (A)
endomat = sparse(total_num_rows ,nrh +nre);

%% Filling the combined stochiometric matrix (A)
%   HOST          GUEST
% [ S_ext ] | [     0     ] [>] [(Host's rhs_ext_lb ) + (Guest's rhs_ext_lb )]-> REGION b1
% [ S_ext ] | [     0     ] [<] [(Host's rhs_ext_ub ) + (Guest's rhs_ext_ub )]-> REGION b2  
% [   0   ] | [S_unmapped ] [=] [0]                                           -> REGION b3  
% [ S_int ] | [ S_ext2int ] [=] [0]                                           -> REGION b4 
% [   0   ] | [   S_int   ] [=] [0]                                           -> REGION b5  
% [   1   ] | [    -1     ] [=] [0]   -> REGION for s.t (Host's growth) = (Guest's growth)

% REGION b1: e compartment lower bound
endomat(b1 ,1:nrh) = hemodel1.S_ext;
% endomat(b1 ,nrh+1:nre+nrh)= 0;        

% REGION b2: e compartment upper bound
endomat(b2 ,1:nrh) = hemodel1.S_ext;
%endomat(b2 ,nrh+1:nre+nrh) = 0;        

% REGION b3:
%endomat(b3 ,1:nrh) = 0;                
endomat(b3 ,nrh+1:nre+nrh) = hemodel2.S_unmapped;

% REGION b4: host's intra. and guest's extrac. region 
endomat(b4 ,1:nrh) = hemodel1.S_int;
endomat(b4 ,nrh+1:nre+nrh) = hemodel2.S_ext2int;

% REGION b5: c compartment endo
%endomat(b5 ,1:nrh) = 0;                
endomat(b5 ,nrh+1:nre+nrh) = hemodel2.S_int;






% REGION for optimization
endomat(end ,bmih) = 1;
endomat(end ,bmie+nrh) = -1; % ensures they grow at same rate

% The combined stochiometric matrix (A) must be sparse to run in gurobi
endomodel.A = sparse(endomat);


%% Create .rhs field to use as condition for the combined stochiometric matrix (A) 
% Pre-empty field .rhs
endomodel.rhs = zeros(total_num_rows ,1);

% filing the field .rhs
endomodel.rhs(b1) = hemodel1.rhs_ext_lb + hemodel2.rhs_ext_lb;
endomodel.rhs(b2) = hemodel1.rhs_ext_ub + hemodel2.rhs_ext_ub;
%endomodel.rhs(b3) = zeros(nu, 1); %
%endomodel.rhs(b4) = zeros(ni, 1); %this makes the (host's intrac rhs) + (guest's extrac rhs) = 0
%endomodel.rhs(b5) = zeros(ni, 1); %this makes the (host's intrac rhs) + (guest's intrac rhs) = 0
% endmodel.rhs(end)=0;

%% Create .sense field which defines the condition btw the matrix A and .rhs
endomodel.sense = [repmat('>' ,ne ,1); ...
                   repmat('<' ,ne ,1); ...
                   repmat('=' ,nu ,1);...
                   repmat('=' ,ni ,1);...
                   repmat('=' ,ni ,1);...
                   repmat('=' ,1  ,1)];

%% Update fields for fluxes .lb & .ub
endomodel.lb=[hemodel1.lb;hemodel2.lb];
endomodel.ub=[hemodel1.ub;hemodel2.ub];

% some models have an upper limit on biomass
endomodel.lb(bmih)=0;
endomodel.ub(bmih)=1000;
endomodel.lb(bmie+nrh)=0;
endomodel.ub(bmie+nrh)=1000;
end
