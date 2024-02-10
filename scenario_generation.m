function generated_scenarios=scenario_generation(historical_data,N_scenarios)
            % DoCopula Estimates a copula pdf from the input grouped data and samples from this Nsamples random scenarios
            % copulMethod ([1]: T-copula| [2]: Gaussian)
            % Example: z4=test1.DoCopula(10,2);
            % (test1 is the object)
            if size(historical_data,2) == 1 || size(historical_data,1) == 1
                disp('You gave ungrouped data and this method requires grouped data');
                generated_scenarios = historical_data;
            else
                % ---Estimate pdf of my multivariate input grouped data
                
                rand_var = historical_data';
                    % Transform the data to the copula scale (unit square) using a kernel estimator of the cumulative distribution function.
                    uni_rand_var = zeros(size(historical_data,2),size(historical_data,1));
%                     if objData.figControl == 1
%                         figure; hold on;
                        for i=1:size(historical_data,1)
                            uni_rand_var(:,i)=ksdensity(rand_var(:,i),rand_var(:,i),'function','cdf','Bandwidth',2);
                          %  [f,xi] = ksdensity(rand_var(:,i),'function','pdf','Bandwidth',2);
                       end    
                
                    %-----GAUSSIAN COPULA-----
                    Rho_2 = copulafit('Gaussian',uni_rand_var);
                    copul_rand = copularnd('Gaussian',Rho_2,N_scenarios);
                    % Transform the random sample back to the original scale of the data.
                    for i=1:size(historical_data,1)
                        gen_sce(:,i)=ksdensity(rand_var(:,i),copul_rand(:,i),'function','icdf','Bandwidth',2);
                    end
                    generated_scenarios = gen_sce';
                   
                            end
            end
