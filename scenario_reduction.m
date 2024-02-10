function [load_sel_scen, wind_sel_scen, solar_sel_scen, probVeci]=scenario_reduction(...
 loadscenrank,windscenrank,solarscenrank,loadi,windi,solari,gen_scen_load,gen_scen_wind,...
 gen_scen_solar,cluster_num)
for i=1:1000
    X(i,1)=loadscenrank(loadi(i));
    X(i,2)=windscenrank(windi(i));
    X(i,3)=solarscenrank(solari(i));
end
[i_dx,Centroids] = kmeans(X,cluster_num);

 closest_index = zeros(max(i_dx),1);
            size_iCluster=zeros(max(i_dx),1);
            for iCluster = 1:max(i_dx)
                %# find the points that are part of the current cluster
                currentPointIdx = find(i_dx==iCluster);
                %# find the index (among points in the cluster)
                %# of the point that has the smallest Euclidean distance from the centroid
                %# bsxfun subtracts coordinates, then you sum the squares of
                %# the distance vectors, then you take the minimum
                [~,minIdx] = min(sum(bsxfun(@minus,X(currentPointIdx,:),Centroids(iCluster,:)).^2,2));
                %# store the index into X (among all the points)
                closest_index(iCluster) = currentPointIdx(minIdx);
                size_iCluster(iCluster,1) = length(X(i_dx==iCluster,1));
            end
            
            
            
            
           
            closestPoints = closest_index;
            probVeci =size_iCluster./(size(gen_scen_load,2));
           
            for iCluster = 1:max(i_dx)
                    load_sel_scen(:,iCluster) =  gen_scen_load(:,loadi(closestPoints(iCluster)));
                    wind_sel_scen(:,iCluster) =  gen_scen_wind(:,windi(closestPoints(iCluster)));
                    solar_sel_scen(:,iCluster) =  gen_scen_solar(:,solari(closestPoints(iCluster)));
%                     hold on;
                end
%               