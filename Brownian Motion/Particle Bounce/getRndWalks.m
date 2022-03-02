function coords = getRndWalks(D,cellRadius,nSpins,nTimeSteps,dt)
  %% allocate memory
  
  coords = zeros(3,nSpins,nTimeSteps);
  
  %% usein_boundsll constrants
  DConst = sqrt(2*D*dt);  %tuned parameter to obtain the correct mult by 2 to make 
  radius_mag = cellRadius^2;

  coords(:,:,1) = zeros(3,nSpins);
 

    for ii = 2:nTimeSteps   

        % run steps

        coords(:,:,ii) = coords(:,:,ii-1) + randn(3,nSpins)*DConst;
        
        coords_mag = squeeze(dot(coords(:,:,ii),coords(:,:,ii),1));
        in_bounds  = sqrt(coords_mag)>cellRadius;
        num_out_bounds = sum(in_bounds(:));
        
        if num_out_bounds > 0
            collision_idx = find(in_bounds>0);     

            for collision = 1:length(collision_idx)
                
                trajectory = coords(:,collision_idx(collision),ii)-coords(:,collision_idx(collision),ii-1);
                point_of_intersection = vectorSphereIntersection(cellRadius, coords(:,collision_idx(collision),ii-1),trajectory);
                
                overshoot = coords(:,collision_idx(collision),ii)-point_of_intersection;
                perc_moved = dot(trajectory,overshoot);
                
%                 perc_moved = dot(overshoot,overshoot)/dot(trajectory,trajectory);

                
                new_coord = point_of_intersection - perc_moved*trajectory;

%                 new_coord = point_of_intersection - perc_moved*trajectory;
                coords(:,collision_idx(collision),ii) = reflectedPoint(point_of_intersection, new_coord);

                if dot( coords(:,collision_idx(collision),ii), coords(:,collision_idx(collision),ii)) > (cellRadius^2)
                    fprintf("%d %d\n\n", dot( coords(:,collision_idx(collision),ii), coords(:,collision_idx(collision),ii)),perc_moved);
                    coords(:,collision_idx(collision),ii) = coords(:,collision_idx(collision),ii)*0.999;
                end
            end
            
        end
        
        
        
    end

    
%  Plot of final distribution %%%
    figure(99);hold off 

    for particle = 1:nSpins
        x_ = squeeze(coords(1,particle,:));
        y_ = squeeze(coords(2,particle,:));
        z_ = squeeze(coords(3,particle,:));

        plot3(x_, y_, z_, 'Color', rand(1,3), 'MarkerSize', 9);
        xlabel('x'), ylabel('y'),zlabel('z');
        hold on; 
    end
    hold off
    
end
