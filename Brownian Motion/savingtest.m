% 
% %need to let particle out if it is within the cut out region
% function result = CheckCutOut(origin,x,y,z, sphere_X,sphere_Y,sphere_Z,cutout_disk,cutout_idx)
% %just assume start at origin
%     result = 0;
% 
%     %try testing each of the x,y,x and x,y,z sphere values in a loop
%     for i = 1:size(cutout_idx,1)
%         
%         S_index = cutout_idx(i,find(cutout_idx(i,:)~=0));
%         new_S_index = 1:length(find(squeeze(cutout_disk(i,:,:))==1));
%         
%         sphere_X_val(i,new_S_index) = sphere_X(S_index);
%         sphere_Y_val(i,new_S_index) = sphere_Y(S_index);
%         sphere_Z_val(i,new_S_index) = sphere_Z(S_index);
%         
%     end
%     
%     
%     for num_cutouts = 1:size(cutout_disk,1)
%         
%         for i = 1:length(find(sphere_Z_val(num_cutouts,:)~=0))
%             
%             Sx = sphere_X_val(num_cutouts,i);
%             Sy = sphere_Y_val(num_cutouts,i);
%             Sz = sphere_Y_val(num_cutouts,i);
%             
%             %check if it is meant to be negative or positive in the axis limits
%             if ( (Sx >= 0) && (Sy >= 0) && (Sz >= 0) )
%                 
%                 if ( (Sx <= x) && (Sy <= y) && (Sz <= z) )
%                     result = 1;
%                 end
%                 
%             elseif ( (Sx >= 0) && (Sy >= 0) && (Sz < 0) )
%                 
%                 if ( (Sx <= x) && (Sy <= y) && (Sz >= z) )
%                     result = 1;
%                 end
% 
%             elseif ( (Sx >= 0) && (Sy < 0) && (Sz >= 0) )
%                 
%                 if ( (Sx <= x) && (Sy >= y) && (Sz <= z) )
%                     result = 1;
%                 end
%                 
%             elseif ( (Sx < 0) && (Sy>= 0) && (Sz >= 0) )
%                 
%                 if ( (Sx >= x) && (Sy <= y) && (Sz <= z) )
%                     result = 1;
%                 end      
%                
%             elseif ( (Sx >= 0) && (Sy < 0) && (Sz < 0) )
%                 
%                 if ( (Sx <= x) && (Sy >= y) && (Sz >= z) )
%                     result = 1;
%                 end      
%                 
%             
%             elseif ( (Sx < 0) && (Sy < 0) && (Sz < 0) )
%                 
%                 if ( (Sx >= x) && (Sy >= y) && (Sz >= z) )
%                     result = 1;
%                 end  
%                 
%             end %end of if,elseif,..... statement
%         end %end of all index checks
%     end %end of num_counts 
% 
% end


% %need to let particle out if it is within the cut out region
% function result = CheckCutOut(origin,x,y,z, sphere_X,sphere_Y,sphere_Z,cutout_disk,cutout_idx)
% %just assume start at origin
%     result = 0;
% 
%     %try testing each of the x,y,x and x,y,z sphere values in a loop
%     for i = 1:size(cutout_idx,1)
%         sphere_X_val(i,1:length(find(squeeze(cutout_disk(i,:,:))==1))) = sphere_X(cutout_idx(i,find(cutout_idx(i,:)~=0)));
%         sphere_Y_val(i,1:length(find(squeeze(cutout_disk(i,:,:))==1))) = sphere_Y(cutout_idx(i,find(cutout_idx(i,:)~=0)));
%         sphere_Z_val(i,1:length(find(squeeze(cutout_disk(i,:,:))==1))) = sphere_Z(cutout_idx(i,find(cutout_idx(i,:)~=0)));
%     end
%     
% %     sphere_X = sphere_X(cutout_idx);
% %     sphere_Y = sphere_Y(cutout_idx);
% %     sphere_Z = sphere_Z(cutout_idx);
%     
%     for num_cutouts = 1:size(cutout_disk,1)
%         
%         for i = 1:length(find(sphere_Z_val(num_cutouts,:)~=0))
%             Sx = sphere_X_val(num_cutouts,i);
%             Sy = sphere_Y_val(num_cutouts,i);
%             Sz = sphere_Y_val(num_cutouts,i);
%             
%             %check if it is meant to be negative or positive
%             if ( (Sx >= 0) && (Sy >= 0) && (Sz >= 0) )
%                 
%                 if ( (Sx <= x) && (Sy <= y) && (Sz <= z) )
%                     result = 1;
%                 end
%                 
%                 
%             elseif ( (Sx >= 0) && (sphere_Y_val(num_cutouts,i) >= 0) && (sphere_Z_val(num_cutouts,i) < 0) )
%                 
%                 if ( (sphere_X_val(num_cutouts,i) <= x) && (sphere_Y_val(num_cutouts,i) <= y) && (sphere_Z_val(num_cutouts,i) >= z) )
%                     result = 1;
%                 end
% 
%             elseif ( (sphere_X_val(num_cutouts,i) >= 0) && (sphere_Y_val(num_cutouts,i) < 0) && (sphere_Z_val(num_cutouts,i) >= 0) )
%                 
%                 if ( (sphere_X_val(num_cutouts,i) <= x) && (sphere_Y_val(num_cutouts,i) >= y) && (sphere_Z_val(num_cutouts,i) <= z) )
%                     result = 1;
%                 end
%                 
%             elseif ( (sphere_X_val(num_cutouts,i) < 0) && (sphere_Y_val(num_cutouts,i) >= 0) && (sphere_Z_val(num_cutouts,i) >= 0) )
%                 
%                 if ( (sphere_X_val(num_cutouts,i) >= x) && (sphere_Y_val(num_cutouts,i) <= y) && (sphere_Z_val(num_cutouts,i) <= z) )
%                     result = 1;
%                 end      
%                
%             elseif ( (sphere_X_val(num_cutouts,i) >= 0) && (sphere_Y_val(num_cutouts,i) < 0) && (sphere_Z_val(num_cutouts,i) < 0) )
%                 
%                 if ( (sphere_X_val(num_cutouts,i) <= x) && (sphere_Y_val(num_cutouts,i) >= y) && (sphere_Z_val(num_cutouts,i) >= z) )
%                     result = 1;
%                 end      
%                 
%             
%             elseif ( (sphere_X_val(num_cutouts,i) < 0) && (sphere_Y_val(num_cutouts,i) < 0) && (sphere_Z_val(num_cutouts,i) < 0) )
%                 
%                 if ( (sphere_X_val(num_cutouts,i) >= x) && (sphere_Y_val(num_cutouts,i) >= y) && (sphere_Z_val(num_cutouts,i) >= z) )
%                     result = 1;
%                 end      
%                 
%             else
%                 result = 0;
%                 
%             end %end of if,elseif,..... statement
%         end %end of all index checks
%     end %end of num_counts 
% 
% end