% 2023.06.12 - Lucas S. Souza
% This function 'run_gh_model_uei_Ng(newmodel ,Ng)' 
% is a modified version of 'run_gh_model_ext_int(newmodel)'
% which is a modified Eric's 'runehmodel.m' function
%
% Aim: to incorporate separation of metabolites into ext and int env
%      which inpacts the structure of the Stoichiom matrix .A, .rhs, .sense
% Output: growth rate ('result1') and the modified model ('tempmodel')
%    .A    .sense      .rhs
% [ Ng*S_ext ] [ > ] [ .rhs_ext_lb ] 
% [ Ng*S_ext ] [ < ] [ .rhs_ext_ub ]
% [ Ng*S_int ] [ > ] [ .rhs_int_lb ]
% [ Ng*S_int ] [ < ] [ .rhs_int_ub ]
 

function [result1,tempmodel] = run_gh_model_uei_Ng(newmodel ,Ng)

params = struct();
params.OutputFlag = 0;

%% Defining matrix .A
tempmodel = newmodel;

% tempmodel.A = [newmodel.S_ext; ...
%                newmodel.S_ext; ...
%                newmodel.S_int; ...
%                newmodel.S_int  ...
%                ];


A = [newmodel.S_ext; ...
    newmodel.S_ext; ...
    newmodel.S_int; ...
    newmodel.S_int  ...
    ];
 tempmodel.A = sparse (Ng*A); %sparse (tempmodel.A);


%% Fluxes' low and upper bound (.lb and .ub)
tempmodel.lb = full(newmodel.lb);
tempmodel.ub = full(newmodel.ub);

%% .obj
tempmodel.obj = zeros(size(tempmodel.A,2),1);
tempmodel.obj(newmodel.bmi) = -1;
%% .rhs
tempmodel.rhs = full([newmodel.rhs_ext_lb; ...
                      newmodel.rhs_ext_ub; ...
                      newmodel.rhs_int_lb; ...
                      newmodel.rhs_int_ub; ...
                      ]);
%% .sense
ne = size(newmodel.S_ext,1);
ni = size(newmodel.S_int,1);

tempmodel.sense = [repmat('>',ne,1); ...
                   repmat('<',ne,1); ...
                   repmat('>',ni,1); ...
                   repmat('<',ni,1) ...
                   ];

%% Output
result1 = gurobi(tempmodel,params);

end