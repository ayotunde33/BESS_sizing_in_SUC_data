function grouped_data = grouping_function(ungrouped_data,group_size)
  if size(ungrouped_data,1) >1 && size(ungrouped_data,2) >1
                disp('You gave already ungrouped data');
                grouped_data = ungrouped_data;
  else
            % GroupSamplesBy Goups vector data by a specified group length
            % Example: test1.iniVec=test1.GroupSamplesBy(96);
            % (test1 is the object)

            L=length(ungrouped_data);
            grouped_data=zeros(group_size,L/group_size);
            i_group=1;
            for i=1:group_size:L
                i_groupj =1;
                for j=i:1:i+group_size-1
                    grouped_data(i_groupj,i_group)=ungrouped_data(j);
                    i_groupj = i_groupj+1;
                end
                i_group = i_group+1;
            end
        end
end 

