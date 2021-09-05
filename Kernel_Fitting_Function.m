function LL_k = Kernel_Fitting_Function(k,Obs,Dist,Area,Proportion_Adult_Sampled,Kernel,POWER)


% Grab some parameters from the inputs
Num_Source_Reefs = size(Obs,1)-1;
Num_Sampled_Source_Reefs = sum(Proportion_Adult_Sampled > 0);
Num_Dest_Reefs = size(Obs,2);

% For each reef, create the expected proportions between every reef pair
Proportions = Kernel(Dist, k, POWER);

% Inflate the values by the total area of the source reef
% (larger reefs send out more larvae, all else being equal)
Settlers = Proportions.*repmat(Area,1,size(Proportions,2));

%% ==== CREATED A LIKELIHOOD MATRIX BASED ON THE KERNEL ====

KernelPredictedSettlers = zeros(size(Obs));

% Go through all the source reefs
for i = 1:Num_Source_Reefs

   % What proportion of this source reef (i) was sampled?
   This_PAS = Proportion_Adult_Sampled(i);
   
   % Go through all the destination reefs
   for j = 1:Num_Dest_Reefs

      SettlersFromAssignedReefs = Settlers(i,j);

      % Not all settlers from assigned reefs will be assigned, because not all adults were sampled
      KernelPredictedSettlers(i,j) = SettlersFromAssignedReefs.*(This_PAS.^2 + 2.*This_PAS.*(1 - This_PAS));
      KernelPredictedSettlers(Num_Source_Reefs+1,j) = KernelPredictedSettlers(Num_Source_Reefs+1,j) ...
                                                      + SettlersFromAssignedReefs.*(1-This_PAS).^2;
   end
end

% Loglikelihoods can't handle zeros. Make them very small instead.
KernelPredictedSettlers(KernelPredictedSettlers==0) = 1e-12;

% Normalise values into multinomial probabilities
PredictedProportions = KernelPredictedSettlers./repmat(sum(KernelPredictedSettlers),Num_Source_Reefs+1,1);

LL_k = 0; % Initialise
for i = 1:Num_Dest_Reefs
   % Go through reefs one-by-one: What is the LL of each reef's observations?
   
   ObsVector = Obs(:,i);
   ProbVector = PredictedProportions(:,i);
   
   % Calculate the log likelihood for this particular reef, given the sample that's been observed
   LL_k = LL_k + sum(ObsVector.*log(ProbVector));
   
end

% Invert the likelihood if using a matlab routine that minimises (e.g., FMINCON, FMINSEARCH, FMINBND)
LL_k = -LL_k;




