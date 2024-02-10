function ranked_Scenarios = kantorovich_ranking(generated_scenarios,Nsce)
            % ReduceScenariosBy Create the reduced set of scenarios (either usual or unusual scenarios)
            % UglyOrUsual ([0]: Keep usual scenarios | [1]: Keep unusual scenarios)
            % Example 1: [z1,z2,z3]=test1.ReduceScenariosBy(0,3);
            % (Keep the 3 most usual scenarios = Reject the most 365-3 unusual scenarios)
            % Example 2: [z1,z2,z3]=test1.ReduceScenariosBy(1,365-3); 
            % (Keep the 3 most unusual scenarios)
            % (test1 is the object)
            
                data_sce=generated_scenarios;
            N_sel_scenarios = Nsce;
            
            %----- ITER=1 (find the most similar scenario) ----------------------------
            scen_all=data_sce; % ALL POSSIBLE SCENARIOS
            scen_current=scen_all; % CURRENT POOL OF SCENARIOS
            indexes=1:size(scen_all,2);
            indexes_vector=indexes'; % VECTOR TO KEEP INDEXES
            int_var1=indexes;
            iter_i=1; % ITERATION INDEX
            
            %-----C Matrix creation-----
            C_matrix=zeros(size(scen_all,2),size(scen_all,2)); % COST MATRIX INITIALIZATION
            cost_M=zeros(size(scen_all,2),size(scen_all,2));
            for t=1:size(scen_all,1)
                for s1=1:size(scen_all,2)
                    i=s1;
                    for s2=1:size(scen_all,2)
                        cost_M(i,s2)=abs(scen_current(t,i)-scen_current(t,s2)); % COST UPDATES
                    end
                end
                C_matrix=C_matrix+cost_M;
            end
            %-----Pr vector creation-----
            C_ini=C_matrix;
            Prob(1:size(scen_all,2),1)=1/size(scen_all,2); % SCENARIOS ARE EQUIPROPABLE
            Prob_iter=Prob; % PROBABILITY VECTOR FOR EACH SCENARIO
            dist=C_matrix*Prob_iter; % DISTANCE METRIC VECTOR CALCULATION
            index_min=find(dist==min(dist)); % INDEX OF SCENARIO WITH MIN DISTANCE METRIC
            if size(index_min,1)>1
                index_min=index_min(1); % CHECK FOR DRAWS AND CHOOSE THE 1ST
            end
            Prob_iter(index_min)=0; % PROBABILITY OF THE SCENARIO WITH MIN DISTANCE METRIC
            
            %-----Selected and reJected  sets creation-----
            %             Ws = zeros(size(scen1,1),DesiredScenNmr);
            Ss=scen_current(:,index_min);  % SET OF SELECTED SCENARIOS (Selected)
            Ss_index(1,iter_i)=index_min; % SET OF INDEXES OF SELECTED SCENARIOS
            Ss_index_vec=Ss_index'; % VECTOR OF OF INDEXES OF SELECTED SCENARIOS
            scen_current(:,index_min)=[]; % SELECTED SCENARIO REMOVAL FROM THE POOL
            Rs=scen_current; % SET OF REJECTED SCENARIOS (reJected)
            Rs_index=setdiff(indexes_vector,Ss_index_vec); % SET OF INDEXES OF REJECTED SCENARIOS
            cost_M_new=zeros(size(scen_all,2),size(scen_all,2)); % UPDATE THE COST MATRIX
            %----- END ITER=1 ---------------------------------------------------------
            
            %----- ITERATIONS UPDATES -------------------------------------------------
            
            while size(Ss,2)<N_sel_scenarios
                iter_i=iter_i+1;
                for s2=1:size(scen_all,2)
                    cost_M_new(:,s2)=min(C_matrix(:,int_var1(index_min)),C_matrix(:,s2));
                end
                cost_M_new(int_var1(index_min),:)=C_matrix(int_var1(index_min),:);
                C_matrix=cost_M_new;
                dist= zeros(size(Rs_index,1),1);
                %                 d=[]; % RE-INITIALIZE DISTANCE METRIC VECTOR
                
                int_var1(index_min)=[]; % INTERMEDIATE VARIABLE OF SCENARIO INDEXES HAVING REMOVED THE INDEX OF ITER-1
                for i=1:size(Rs_index,1)
                    dist(i,1)=C_matrix(Rs_index(i),:)*Prob_iter; % RE-CALCUATE THE DISTANCE VECTOR FOR THE REMAINING SCENARIOS
                end
                index_min=find(dist==min(dist)); % INDEX OF NEXT SCENARIO WITH MIN DISTANCE METRIC
                if size(index_min,1)>1
                    index_min=index_min(1); % CHECK FOR DRAWS AND CHOOSE THE 1ST
                end
                Ss=[Ss scen_all(:,int_var1(index_min))]; % INCREASE THE SET OF SELECTED SCENARIOS
                size(Ss,2);
                scen_current(:,index_min)=[]; % SELECTED SCENARIO REMOVAL FROM THE POOL
                Rs=scen_current; % UPDATED SET OF REJECTED SCENARIOS
                Ss_index(1,iter_i)=int_var1(index_min); % UPDATED SET OF INDEXES OF SELECTED SCENARIOS
                Ss_index_vec=Ss_index'; % UPDATED VECTOR OF OF INDEXES OF SELECTED SCENARIOS
                Rs_index=setdiff(indexes_vector,Ss_index_vec); % UPDATED SET OF INDEXES OF REJECTED SCENARIOS
                %Pr_iter(id_min)=0; % SELECTED (REMOVED) SCENARIO HAS ZERO PROBABILITY NOW
                Prob_iter(int_var1(index_min))=0; % SELECTED (REMOVED) SCENARIO HAS ZERO PROBABILITY NOW
            end
            %----- END ITERATIONS UPDATES ---------------------------------------------
            
            %----- RE-ASSIGN PROPABILITIES TO THE SELECTED SCENARIOS ------------------
            % Initialize the prob of each selected scenario as the default prob
            prob_s = zeros(size(Ss,2),1);
            for s=1:size(Ss,2)
                selected_scen=Ss_index(s);
                prob_s(s)=Prob(selected_scen);
            end
            % Transfer prob of the rejected scenarios to the closest selected ones
            distance = zeros(size(Ss,2),1);
            for j=1:size(Rs,2)
                for s=1:size(Ss,2)
                    selected_scen=Ss_index(s);
                    distance(s)=C_ini(Rs_index(j),selected_scen);
                end
                mindist_s=find(distance==min(distance));
                prob_s(mindist_s)=prob_s(mindist_s)+Prob(Rs_index(j));
            end
            %----- END RE-ASSIGN PROPABILITIES ----------------------------------------
            
            
            %------- MAIN OUTPUT ------------------------------------------------------
            % test1.iniVec = UnGroupSamples(test1);  (Ungroup to initial vector)
            %{
            Wj_obj = DataX(Rs);
            Ws_obj = DataX(Ss);
            if UglyOrUsual == 0
                % OUTPUT: SELECTED SCENARIOS (The most usual-similar ones)
                %                 selScen=days2ScenVec(Ws);
                selScen = UnGroupSamples(Ws_obj);
                reassignedProbs=prob_s;
                orderedIndexes=Ss_index;
            else
                % OUTPUT: REJECTED SCENARIOS (The most unusual-different ones)
                %                 selScen=days2ScenVec(Wj);
                selScen = UnGroupSamples(Wj_obj);
%                 prob(1:size(Ws,2))=1/size(Wj,2);
                reassignedProbs(1:size(Rs,2),1)=1/size(Rs,2);
                orderedIndexes=Rs_index;
                %------- END OUTPUT -------------------------------------------------------
            end
            %}
        %    reassignedProbs=prob_s;
                ordered_Indexes=Ss_index;
        
%        
        
            % RankScenarios Rank scenarios of the dataset based on the Kantorovich similarity index
            % Example: k = RankScenarios(test1);
            % (test1 is the object)
            
%             [~,~,orderedIndexes]=objData.ReduceScenariosBy(0,size(objData.iniVec,2)-1);
            
%             [~,~,orderedIndexes]=objData.ReduceScenariosBy(0,364);

            score_scale = 1;
            i_scen=1;
            Linear_score = linspace(score_scale,0,size(ordered_Indexes,2));
            ranked_Scenarios = zeros(size(Linear_score,2),1);
            while i_scen<=size(Linear_score,2)
                for i=1:size(Linear_score,2)
                    if ordered_Indexes(i) == i_scen
                        ranked_Scenarios(i_scen) = Linear_score(i);
                    end
                end
                i_scen=i_scen+1;
            end
        end    
        
%}