function [artifact_components, r_squared, beta_array, residuals_array] = TMSfMRIart(x_mat_filename,y_mat_filename, threshold, model_predictors)
% TMSfMRI_explainVariance 
% > uses timing & intensity of TMS to flag specific MELODIC components for removal via a simple multiple regression model
% > JMS 2017
%
% >>>>> INPUT VALUES <<<<<
%
% > x_mat_filename: file containing variable named 'x_mat'
%   > x_mat is a matrix of predictors, e.g. Constant, StimA_time, StimB_time, StimA_time+Lag1, StimB_time+Lag1 
%   > n-by-m: n timepoints, m predictors
%   > the first column should be a vectors of 1's to allow for an intercept in regression model
%   > the other columns are generally 0's and 1's, but don't have to be
%   > The sample matrices includes these predictors:
%       > Column 1: vector of 1's for intercept
%       > Column 2: 1 for TRs that contain TMS at intensity 1, 0 for all other TRs
%       > Column 3: 1 for TRs that contain TMS at intensity 2, 0 for all other TRs
%       > Column 4: same as column 2 but shifted forward one TR in time (so now the TR following TMS is identified) [Lag1]
%       > Column 5: same as column 3 but shifted forward one TR in time (so now the TR following TMS is identified) [Lag1]
%       > Column 6: same as column 2 but shifted forward two TRs in time [Lag2]
%       > Column 7: same as column 3 but shifted forward two TRs in time [Lag2]
%       > Column 8: same as column 2 but shifted *backward* one TR in time [Lag-1 --- can be used as a control]
%       > Column 9: same as column 3 but shifted *backward* one TR in time [Lag-1 --- can be used as a control]
% 
% > y_mat_filename: file containing variable named 'y_mat'  
%   > y_mat is a matrix of components from FSL's MELODIC 
%   > n-by-m: n timepoints, m components
%
% > threshold: threshold for variance explained (R^2) to flag a component as an artifact; in the range 0 to 1
% 
% > model_predictors: columns of x_mat to use as predictors
%   > should always include Column 1 (should be a vector of 1's) of X for the intercept
%   > pick other columns to include as predictors in regression analysis, e.g. [1 2 3 4 5] or [1 4 5]
%
% > NOTES:
%   > this does not currenly use spatial information to flag components for removal (over and above the inherent spatiotemporal nature of MELODIC)
%
%%%%% sample command: 
% [artifact_components, r_squared] = TMSfMRIart('run1x','run1y_P003', .4, [1 2 3 4 5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load variables
load(x_mat_filename); load(y_mat_filename);

% check that X and Y contain an equal # of timepoints before trying regression
if size(x_mat,1) ~= size(y_mat,1)
    error('ERROR: X and Y do not contain equal number of timepoints');
end


% Multiple Linear Regression: Loop through y_mat column by column to find model fits for each MELODIC component
for loop1=1:size(y_mat,2)
    [beta,~,residuals,~,stats] = regress(y_mat(:,loop1),x_mat(:,model_predictors));        
    r_squared(loop1)=stats(1); % the r-squared is the first element of 'stats', which also incldes F, p, etc
    beta_array(:,loop1)=beta;
    residuals_array(:,loop1)=residuals;
end


% print to command window the current configuration and the components flagged as artifacts
disp(' ');
disp(['Number of model parameters (incl intercept) is: ' num2str(size(model_predictors,2))]);
disp(['Maximum model fit (R^2) across all components: ' num2str(round(max(r_squared)*100)) '%']);
disp(['Average model fit (R^2) across all components: ' num2str(round(mean(r_squared)*100)) '%']);
disp(['Threshold: ' num2str(threshold*100) '% of variance']);
disp(' ');
artifact_components=[];
for loop2=1:size(r_squared,2)
    if r_squared(loop2) >= threshold
        disp(['Likely artifact: component #' num2str(loop2) ' (' num2str(round(r_squared(loop2)*100)) '%)'])
        artifact_components=[artifact_components loop2];
    end
end
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potentially useful figures to generate at the end of the script:

% >>> Plot all of the r_squared values; useful to see spread and ensure threshold is appropriate
if 1 % set to 1 to enable
    figure; EVbar=bar(r_squared);
    title('Explained Variance for Each Component');
    xlabel('Components');
    ylabel('Explained Variance (R^2)');
end

% >>> Plot the betas for each selected component
if ~isempty(artifact_components) && 1 % set to 1 to enable
    figure; bar(transpose(beta_array(:,artifact_components)));
    title('Regression Coefficients for Each Artifact Component');
    xlabel('Artifact Components');
    ylabel('Beta Weight');
end

% >>> Generate the best-fit regression line for any component; plot that line, the corresponding residual, and the raw component timecourse:
%   > This can be useful for exploring how the residual structure of a given component changes depending on which regressors are included.
if ~isempty(artifact_components) && 0 % set to 1 to enable
    % Pick an artifact component for which to do this:
    componentNum=artifact_components(3); 
    % Or, pick an arbitrary component:
    % componentNum=1; 
    figure;
    regressLine=x_mat(:,model_predictors)*beta_array(:,componentNum);
    hold on
    plot(regressLine,'color',[0 1 0]); % plot regression line in green
    plot(residuals_array(:,componentNum),'color',[1 0 0]); % plot residuals in red
    plot(y_mat(:,componentNum),'color',[0 0 1]); % plot raw component timecourse in blue
    hold off
end


% >>> Put a stop here in Editor to stop function and have access to internal variables
stopVar=1;


