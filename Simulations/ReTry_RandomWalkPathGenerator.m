function ReTry_RandomWalkPathGenerator(nSet,nSpins,Vradius,nTimeSteps,diffusionGradient2_loc,D,n,dt)
    %%%%%%%%%%%%% COORDS SETUP %%%%%%%%%%%%%
    Coords = zeros(length(Vradius),3,nSpins,nTimeSteps); %particle start loc is assume 0,0,0

    for set = 1:nSet %20 sets of 500 for 10000 spins

        for radius_bounds = 1:length(Vradius)
            j = 1;
            fprintf("RADIUS %d ----- SET %d  \n",radius_bounds, set);
            %starting spin locations
            r = StartingPos(Vradius(radius_bounds),nSpins);

            Coords(radius_bounds,1,:, j) = r(1,:)';
            Coords(radius_bounds,2,:, j) = r(2,:)';
            Coords(radius_bounds,3,:, j) = r(3,:)';

                for j = 2:nTimeSteps   

                    %dxdydz step distributions
                    sd__ = 0.58;%sqrt(STOP_TIME/numberOfSteps);
                    rnd_x = sd__.*randn(1,nSpins)*sqrt(2*n*D*dt);
                    rnd_y = sd__.*randn(1,nSpins)*sqrt(2*n*D*dt);
                    rnd_z = sd__.*randn(1,nSpins)*sqrt(2*n*D*dt);


                    temp_r = r + [rnd_x;rnd_y;rnd_z]; %adding the step

                    %%% Bounds checking %%%
                    r = InBounds(Vradius(radius_bounds),temp_r(1,:),temp_r(2,:),temp_r(3,:),rnd_x,rnd_y,rnd_z);

                    Coords(radius_bounds,1,:, j) = r(1,:)';
                    Coords(radius_bounds,2,:, j) = r(2,:)';
                    Coords(radius_bounds,3,:, j) = r(3,:)';

                    % condition is true if it has reached the read point
                    if (j > diffusionGradient2_loc(end))
                        break;
                    end

                end

     %% Plot of final distribution %%%
%             figure; 
%     
%             for particle = 1:nSpins
%                 x_ = squeeze(Coords(radius_bounds,1,particle,1:diffusionGradient2_loc(end)));
%                 y_ = squeeze(Coords(radius_bounds,2,particle,1:diffusionGradient2_loc(end)));
%                 z_ = squeeze(Coords(radius_bounds,3,particle,1:diffusionGradient2_loc(end)));
%                 
%                 plot3(x_, y_, z_, 'Color', rand(1,3), 'MarkerSize', 9);
%                 xlabel('x'), ylabel('y'),zlabel('z');
%                 hold on; 
%             end
%             hold off
%             pause;

        end %end of one radius cycle value
        
        file_path = "C:\Users\s4427550\3D_Coords"; %         save(sprintf('3D_Coords/Coords%d.mat',set),'Coords');

        save(sprintf('%s%d.mat',file_path,set),'Coords');
    end % end of each set of spins
end


function result = InBounds(r,x1,y1,z1,dx,dy,dz)
    %result is dim = 3 x nSpins
    origin = [0 0 0];
    result = [x1;y1;z1];

    mag_ = sqrt((abs(x1)-origin(1)).^2 + (abs(y1)-origin(2)).^2 + (abs(z1)-origin(3)).^2);
    test_res = (mag_>=r);
    temp = max(test_res);
    temp = max(temp);
    
    if temp == 1
        test = 1;
        D = 3e-9;
        n = 3;
        dt = 10.0e-6;
        sd__ = 0.58;%sqrt(STOP_TIME/numberOfSteps);

        index_list_ = find((mag_ >= r) > 0);
        counter = 0;

        while( test == 1 )
            if (counter == 0)
                %dxdydz step distributions
                rnd_x = sd__.*randn(1,length(index_list_))*sqrt(2*n*D*dt);
                rnd_y = sd__.*randn(1,length(index_list_))*sqrt(2*n*D*dt);
                rnd_z = sd__.*randn(1,length(index_list_))*sqrt(2*n*D*dt);

                pos_0 = [(x1(index_list_)-dx(index_list_));...
                         (y1(index_list_)-dy(index_list_));...
                         (z1(index_list_)-dz(index_list_))];
                      
                test_res = pos_0 + [rnd_x;rnd_y;rnd_z];

                mag = sqrt((abs(test_res(1,:))-origin(1)).^2 + (abs(test_res(2,:))-origin(2)).^2 + (abs(test_res(3,:))-origin(3)).^2);
                temp = (mag>=r);
                
                %update the ones that are now correct
                index_list = find((mag < r) > 0); 
                
                if (max(index_list)==1)
                    result(:,index_list_(index_list)) = test_res(:,(index_list));
                    x1(index_list_(index_list)) = result(1,(index_list_(index_list))); 
                    y1(index_list_(index_list)) = result(2,(index_list_(index_list))); 
                    z1(index_list_(index_list)) = result(3,(index_list_(index_list)));
                end
                %test if code should stop
                test = max(temp);
                if (test == 0)
                    break;
                end 
                
                counter = 1;
            else
                
                %find index values that need to be re-generated
                mag = sqrt((abs(result(1,:))-origin(1)).^2 + (abs(result(2,:))-origin(2)).^2 + (abs(result(3,:))-origin(3)).^2);
                index_list_ = find((mag >= r) > 0); 
                 
                %dxdydz step distributions
                rnd_x = sd__.*randn(1,length(index_list_))*sqrt(2*n*D*dt);
                rnd_y = sd__.*randn(1,length(index_list_))*sqrt(2*n*D*dt);
                rnd_z = sd__.*randn(1,length(index_list_))*sqrt(2*n*D*dt);
                
                %remove prpevious delta_pos to add ned one on
                pos_0 = [(x1(index_list_)-dx(index_list_));...
                         (y1(index_list_)-dy(index_list_));...
                         (z1(index_list_)-dz(index_list_))];
                      
                                        
                test_res = pos_0 + [rnd_x;rnd_y;rnd_z];
                
                %test if new step is within the magnitude conditions
                mag = sqrt((abs(test_res(1,:))-origin(1)).^2 + (abs(test_res(2,:))-origin(2)).^2 + (abs(test_res(3,:))-origin(3)).^2);
                temp = (mag>=r);
                
                %update the ones that are now correct
                index_list_update = find((mag < r) > 0); 
                if (max(index_list_update)>0)
                    result(:,index_list_(index_list_update)) = test_res(:,(index_list_update));
%                     x1(index_list_(index_list_update)) = test_res(1,(index_list_update)); 
%                     y1(index_list_(index_list_update)) = test_res(2,(index_list_update)); 
%                     z1(index_list_(index_list_update)) = test_res(3,(index_list_update));
                end
                
                %if no values are above r, then continue code
                test = max(temp);
                if (test == 0)
                    break;
                end 
                   
            end % end of if counter == 1 || counter == 0
        end % end of while
    end % end of (temp == 1)
end


function point_in_sphere = StartingPos(radius,nSpins)
 
	n = nSpins;
    r = (radius-1.0e-9) * (((rand(1,n))').^3);
    costheta = 2*rand(n,1)-1;
    sintheta = sqrt(1-costheta.^2);
    phi=rand(n,1)*(2*pi);
    x=r.*sintheta.*cos(phi);
    y=r.*sintheta.*sin(phi);
    z=r.*costheta;


    point_in_sphere = [x';...
                       y';...
                       z'];
%     figure;               
%     plot3(x,y,z,'.');
%     xlabel('x'), ylabel('y'),zlabel('z');

end


%Draws the plot for the particles trajetories 
function DrawParticleGraph(xCoords,yCoords,zCoords,nSpins,fig_num)
    figure(fig_num)
    for spin = 1:nSpins
        plot3(xCoords(spin,:),yCoords(spin,:),zCoords(spin,:),'Color', rand(1,3), 'MarkerSize', 9);
        hold on; 
    end
    xlabel('x'), ylabel('y'),zlabel('z');
    title('Final location of particle')
    drawnow;
end



