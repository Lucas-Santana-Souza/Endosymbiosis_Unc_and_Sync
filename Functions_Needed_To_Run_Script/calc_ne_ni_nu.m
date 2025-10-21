% 2023.09.19 - Lucas S. Souza
% function to return …
% I create this function so that ne, ni, nu does not need to be calculated for every parfor loop. 
% I can do that because within each database (Agora and CarveMe) the ne, ni, and nu have a same number.

function [ne, ni, nu] = calc_ne_ni_nu (filedir)
    % Load the host model
    % Change 1: Replaced EVAL with LOAD to load the host model, avoiding "transparency violation error"
    host_model_path = fullfile(filedir, ['model1.mat']);
    ehmodel1 = load(host_model_path, 'metabolic_model');
    ehmodel1 = ehmodel1.metabolic_model;


    % Parameters used in the function 'create_holobiont_default_rhs'
    ne = size(ehmodel1.S_ext      ,1); % ne -> # of extracelular mapped metabolites
    ni = size(ehmodel1.S_int      ,1); % ni -> # of intracelullar metabolites 
    nu = size(ehmodel1.S_unmapped ,1); % nu -> # of extracelular unmapped metabolites 

end