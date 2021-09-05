function Fit_Bolbo_Data()

% ================================================
%   Parentage fitting procedures. Contact michael.bode@qut.edu.au with questions, errors, etc.
% ================================================

% % %% ======================== DATA INPUTS ========================

% First, read in the parentage matrix.
load('SolomonsParentageMatrix.mat','Parentage_matrix')
Assignments = Parentage_matrix;

% Second, read in the proportion of adults on each reef that were sampled.
load crcb_domain_Solomons
Adult_sample_proportions = zeros(length(reef_centroid),1);
% We're going to assume that all adults at the adult sample locations were sampled
Adult_sample_proportions(find(sum(Parentage_matrix(1:end-1,:),2)>0)) = 1;

% Third, read in a distance matrix showing the distances between all
load StoredDistanceMatrix DistanceMatrix
Distances = DistanceMatrix;

% Fourth, how large are all the reefs?
Reef_sizes = reef_centroid(:,3);

% Fifth, read in a list of the reefs that were sampled (specifically, for juveniles)
Sampled_reefs_J = find(sum(Parentage_matrix(1:end-1,:))>0);
Sampled_reefs_A = find(sum(Parentage_matrix(1:end-1,:),2)>0);

% Sixth, read in the locations of each of the reefs
Centroids = reef_centroid(:,1:2);

%% CUT OUT THE HARD STUFF
Distances = Distances(:,Sampled_reefs_J);
Assignments = Assignments(:,Sampled_reefs_J);

save FastFittingDataset
clear all
load FastFittingDataset Distances Centroids Adult_sample_proportions Assignments Reef_sizes Sampled_* *_table

Distances = Distances.*110;
Adult_sample_proportions = 0.5.*Adult_sample_proportions;

% Remove the reef area of juvenile habitat, since it can't generate recruits
Reef_sizes(length(A_table)+1:end) = 0;

% Define the generalised Gaussian functions that we're fitting here
F  = @(x,k,theta)    exp(k).*theta.*exp(-(exp(k)*x).^theta)/gamma(1/theta);
FM = @(x,k,theta) x.*exp(k).*theta.*exp(-(exp(k)*x).^theta)/gamma(1/theta);

% Run through the list of potential kernels and fit each one to the data
Theta_list = [1 2 3 0.5];
LowerBound = -14;
UpperBound =  0;
for th = 1:length(Theta_list);
    [Best_k(th),LL_k(th)] = fminbnd(@Kernel_Fitting_Function,LowerBound,UpperBound,[],... % These are the search input parameters
        Assignments,Distances,Reef_sizes,Adult_sample_proportions,F,Theta_list(th)); % These are the extra parameters needed by the function
end

% Identify the best fit from the candidate kernels
[~,Best_kernel] = min(LL_k);
Best_fit_k = Best_k(Best_kernel);

% Create bootstrap confidence bounds around the best estimate by re-sampling at the destination patch scale
Num_Dest_Reefs = length(Sampled_reefs_J);
for b = 1:250
    if mod(b,10) == 0; disp(b); end
    % First resample the assignment data, as well as the proportion of each reef sampled
    Bootstrap_resample = randsample(Num_Dest_Reefs,Num_Dest_Reefs,1);
    Assignments_b = Assignments(:,Bootstrap_resample);
    Distances_b = Distances(:,Bootstrap_resample);
    
    % Generate lognormally distributed abundance estimates
    M = 1;
    SIGMA = 0.5;
    VAR = SIGMA^2;
    mu = log(M^2/sqrt(VAR+M^2));
    sig = sqrt(log(VAR/M^2 + 1));
    Reef_sizes_b = Reef_sizes.*lognrnd(mu,sig,size(Reef_sizes));

    % Re-fit the kernel of the correct shape to this data
    [Bootstrap_k(b),Bootstrap_LL_k(b)] = fminbnd(@Kernel_Fitting_Function,LowerBound,UpperBound,[],... % These are the search input parameters
        Assignments_b,Distances_b,Reef_sizes_b,Adult_sample_proportions,F,Theta_list(Best_kernel)); % These are the extra parameters needed by the function
    
end
Confidence_bounds = quantile(Bootstrap_k,[0.025 0.5 0.975])

% Calculate the mean dispersal distance by integrating the best fit kernel,
% and the confidence intervals from the quantiles of the bootstrap re-sampled fits.
MDD   =  integral(@(x)FM(x,Best_fit_k,Theta_list(Best_kernel)),0,5000);
MedDD   =  integral(@(x)FM(x,Confidence_bounds(2),Theta_list(Best_kernel)),0,5000);
CI_DD = [integral(@(x)FM(x,Confidence_bounds(1),Theta_list(Best_kernel)),0,5000) ...
    integral(@(x)FM(x,Confidence_bounds(3),Theta_list(Best_kernel)),0,5000)];

disp(['Mean dispersal distance is ' num2str(MDD,3) ' km'])
disp(['Median dispersal distance is ' num2str(CI_DD(2),3) ' km'])
disp(['with 95% confidence intervals of [' num2str(CI_DD(1),3) ', ' num2str(CI_DD(2),3) '] km'])

if Best_fit_k == LowerBound
    disp('ERROR: Reduce the lower bound passed to function FMINBND')
elseif Best_fit_k == UpperBound
    disp('ERROR: Increase the upper bound passed to function FMINBND')
end

save OUTPUTS_random_size
Plot_kernel







