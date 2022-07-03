function [AUC, beta] = run_ENGLM(diagnosis_vector,rsFC_predictors)
% run_ENGLM: Generate an elastic-net regularized general linear model for
% predicting diagnostic status from rsFC predictor variables. 

%Inputs: 
%    - diagnosis_vector: single column vector assigning the diagnostic
%    status (1 or 0) for each subject (rows)
%    - rsFC_predictors: subject-by-rsFC feature matrix containing
%    functional connectivity predictors (columns) for each subject (rows)

%Outputs: 
%    - AUC: Area under ROC curve representing GLM performance at different
%    lambda values for bootstrap iteration
%    - beta: GLM beta weights assigned to each predictor at different
%    lambda values for bootstrap iteration

%% 

%define train and test indices
sample_size = 217;
test_sample_size = round(sample_size/3);
train_sample_size = sample_size - test_sample_size;
bootstrap_count = 100;

%define GLM settings 
set_alpha = 0.5;

%pre-allocate 
AUC = cell(1,bootstrap_count);
B = cell(1,bootstrap_count);
Fitinfo = cell(1,bootstrap_count);

%run GLM
for i = 1:bootstrap_count
       
    %define random train and test subsample indices
    train_sample_idx = randsample(size(rsFC_predictors,1), train_sample_size);
    test_sample_idx = randsample(find(~ismember(1:size(rsFC_predictors,1),train_sample_idx)), test_sample_size);
    
    %generate train and test datasets
    train_predictor_var = rsFC_predictors(train_sample_idx,:);
    train_response_var = diagnosis_vector(train_sample_idx);
    test_predictor_var = rsFC_predictors(test_sample_idx,:);
    test_response_var = diagnosis_vector(test_sample_idx);
    
    %run train dataset through GLM
    [temp_B,temp_Fitinfo] = lassoglm(train_predictor_var,train_response_var, 'binomial','Alpha',set_alpha);

    B{i} = temp_B;
    Fitinfo{i} = temp_Fitinfo;    

    temp3 = zeros(1,size(temp_B,2));

    %compute sensitivity and specificity for held-out test sample
    for ii = 1:size(temp_B,2)
        [~,~,~, temp3(ii)] = perfcurve(test_response_var, test_predictor_var * temp_B(:,ii), 1);
    end

    AUC{i} = temp3;
    
end


end

