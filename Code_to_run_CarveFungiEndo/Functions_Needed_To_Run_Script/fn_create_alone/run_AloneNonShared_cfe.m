%% 2026.05.12 - Lucas Santana Souza
%
% Filling the combined stochiometric matrix (H)
%% Structure        
%  H stoichiometric matrix          
% [ S_ext ]  [>] [(Host's rhs_ext_lb) ]-> REGION b1
% [ S_ext ]  [<] [(Host's rhs_ext_ub) ]-> REGION b2  
% [ S_intC ] [=] [0]                   -> REGION b3 
% [ S_intO ] [=] [0]                   -> REGION b4 

function [result1,ConstructedModel] = run_AloneNonShared_cfe(inputModelH ,ne ,nic ,nio ,nu)

params = struct();
params.OutputFlag = 0;

%% Properties/dimentions

nrh = size(inputModelH.lb,1);   % number of columns/reactions of H stoichoimetric matrix
%nre = size(inputModelE.lb,1);   % number of columns/reactions of H stoichoimetric matrix

% find biomass index
bmih = inputModelH.bmi;
%bmie = inputModelE.bmi;

%% construct H Stoichiometric matrix for E and C compartments

%% Defining metabolite regions to be used at the combined stochiometric matrix (A) and combined rhs (.rhs)
b1 = [1:ne];               % b1: row indices corresponding to A for extracellular lb
b2 = [1:ne]  + ne;         % b2: row indices corresponding to A for extracellular ub
b3 = [1:nic] + 2*ne;       % b3: row indices corresponding to A for intracellular  (citoplasm)
b4 = [1:nio] + 2*ne + nic; % b4: row indices corresponding to A for intracellular  (other compartments)

total_num_rows = 2*ne + nic + nio;
total_num_cols = nrh;

%% Pre-empty the combined stochiometric matrix (A)
ConstructedMatrix = sparse(total_num_rows ,total_num_cols);


%% Filling the combined stochiometric matrix (A)
%   HOST          
% [ S_ext ]  [>] [(Host's rhs_ext_lb)]-> REGION b1
% [ S_ext ]  [<] [(Host's rhs_ext_ub)]-> REGION b2  
% [ S_intC ] [=] [0]                  -> REGION b3 
% [ S_intO ] [=] [0]                  -> REGION b4 


% REGION b1: e compartment lower bound
ConstructedMatrix(b1 ,1:nrh) = inputModelH.S_ext;
% endomat(b1 ,nrh+1:nre+nrh)= 0;        

% REGION b2: e compartment upper bound
ConstructedMatrix(b2 ,1:nrh) = inputModelH.S_ext;
%endomat(b2 ,nrh+1:nre+nrh) = 0;         

% REGION b3: host's intra. region C
ConstructedMatrix(b3 ,1:nrh) = inputModelH.S_intC;
%ConstructedMatrix(b3 ,nrh+1:nre+nrh) = inputModelG.S_ext2int;

% REGION b4: host's intra. region O
ConstructedMatrix(b4 ,1:nrh) = inputModelH.S_intO;
%ConstructedMatrix(b4 ,nrh+1:nre+nrh) = 0;

% The combined stochiometric matrix (A) must be sparse to run in gurobi
ConstructedModel.A = sparse(ConstructedMatrix);

%% Create .rhs field to use as condition for the combined stochiometric matrix (A) 

% Pre-empty field .rhs
ConstructedModel.rhs = zeros(total_num_rows ,1);

% filing the field .rhs
 ConstructedModel.rhs(b1) = inputModelH.rhs_ext_lb; %+ inputModelE.rhs_ext_lb;
 ConstructedModel.rhs(b2) = inputModelH.rhs_ext_ub; %+ inputModelE.rhs_ext_ub;
%ConstructedModel.rhs(b3) = zeros(nic, 1); %this is rhs for region S_intC
%ConstructedModel.rhs(b4) = zeros(nio, 1); %this is rhs for region S_intO

%% Create .sense field which defines the condition btw the matrix A and .rhs
ConstructedModel.sense = [repmat('>' ,ne  ,1);... % -> REGION b1
                          repmat('<' ,ne  ,1);... % -> REGION b2
                          repmat('=' ,nic ,1);... % -> REGION b3
                          repmat('=' ,nio ,1)];   % -> REGION b4

%% Fluxes' low and upper bound (.lb and .ub)
ConstructedModel.lb = full(inputModelH.lb);
ConstructedModel.ub = full(inputModelH.ub);

% some models have an upper limit on biomass
ConstructedModel.lb(bmih)=0;
ConstructedModel.ub(bmih)=1000;
%ConstructedModel.lb(bmie+nrh)=0;
%ConstructedModel.ub(bmie+nrh)=1000;

%% Create .obj field, which determines what is maximized

c_obj = zeros(total_num_cols ,1); % It should be of the size of all columns in overal stoic. matrix 
                                  % (here because it is only one species it is just the number of fluxes of
                                  % that species, if there would be two species would be the sum of number of
                                  % fluxes of both species)
c_obj(bmih) = -1;    % because the default in Gurobi is minimize, I have to put '-1'
ConstructedModel.obj = c_obj; 


%% Output
result1 = gurobi(ConstructedModel,params);

end