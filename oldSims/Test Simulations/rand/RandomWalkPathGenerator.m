function RandomWalkPathGenerator(nSet,nSpins,Vradius,nTimeSteps,diffusionGradient2_loc,D,n,dt)
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
            figure; 
    
            for particle = 1:nSpins
                x_ = squeeze(Coords(radius_bounds,1,particle,1:diffusionGradient2_loc(end)));
                y_ = squeeze(Coords(radius_bounds,2,particle,1:diffusionGradient2_loc(end)));
                z_ = squeeze(Coords(radius_bounds,3,particle,1:diffusionGradient2_loc(end)));
                
                plot3(x_, y_, z_, 'Color', rand(1,3), 'MarkerSize', 9);
                xlabel('x'), ylabel('y'),zlabel('z');
                hold on; 
            end
            hold off
            pause;

        end

        save(sprintf('3D_Coords/Coords%d.mat',set),'Coords');
    end
end


function result = InBounds(r,x1,y1,z1,dx,dy,dz)
    %result is dim = 3 x nSpins
    
    origin = [0 0 0];
    result = [x1;y1;z1];

    mag = sqrt((abs(x1)-origin(1)).^2 + (abs(y1)-origin(2)).^2 + (abs(z1)-origin(3)).^2);
    test_axis_x = sqrt((abs(x1)-origin(1)).^2 + (abs(y1-dy)-origin(2)).^2 + (abs(z1-dz)-origin(3)).^2);
    test_axis_y = sqrt((abs(x1-dx)-origin(1)).^2 + (abs(y1)-origin(2)).^2 + (abs(z1-dy)-origin(3)).^2);
    test_axis_z = sqrt((abs(x1-dx)-origin(1)).^2 + ((abs(y1-dy))-origin(2)).^2 + (abs(z1)-origin(3)).^2);
    test_axis = [test_axis_x;test_axis_y;test_axis_z];
    test_res = (test_axis>=r);
    temp = max(test_res);
    temp = max(temp);
    
% 	if (max(test_axis_x >= r))
%         index_list = find((test_axis_x >= r) > 0);
%         x1(index_list) = x1(index_list)-dx(index_list);
%     end
%     
%     if (max(test_axis_y >= r))
%         index_list = find((test_axis_y >= r) > 0);
%         y1(index_list) = y1(index_list)-dy(index_list);
%     end
%     
%     if (max(test_axis_z >= r))
%         index_list = find((test_axis_z >= r) > 0);
%         z1(index_list) = z1(index_list)-dz(index_list);
%     end
    
    
    result = [x1;y1;z1];

    
    
    if (temp == 1)
        
        need_to_mirror = sum(test_res);
        index_list = find(need_to_mirror>0);
    
        mirrored = mirror_trajectory([(x1-dx);(y1-dy);(z1-dz)], [dx;dy;dz],index_list,r);
        result = mirrored;
    end
    
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



%% %%% from other file %%%

% 
% TODO:
% Should probably vectorise this
function result = mirror_trajectory(previous_xyz,current_dxdydz, index_list,radius)
    
    x0 = previous_xyz(1,:);
    y0 = previous_xyz(2,:);
    z0 = previous_xyz(3,:);
    
    x1 = x0 + current_dxdydz(1,:);
    y1 = y0 + current_dxdydz(2,:);
    z1 = z0 + current_dxdydz(3,:);
    
    result = [x1;y1;z1];
    
    vertex_x0y0z0 = [x0',y0',z0']; 
    vertex_x1y1z1 = [x1',y1',z1'];
    
%     figure;
    for i =1:length(index_list)
        
        V = [vertex_x0y0z0(index_list(i),:);...
             vertex_x1y1z1(index_list(i),:)];
         
%         plot3(V(:,1),V(:,2),V(:,3),'k.-','MarkerSize',25,'Color', rand(1,3)); 
%         xlabel('x'), ylabel('y'),zlabel('z');
%         hold on
        
        % NEW PART
        syms r [1 3] 
        f = r*r.';
        
        feqn = (f == radius); %radius is 14^(0.5)
        fgrad = gradient(f,r);

        size(fgrad);

        % Define the equation for the tangent plane. Use the subs function to
        % evaluate the gradient at the point r0
        r0 = vertex_x0y0z0(index_list(i),:);%[0.0911,0.1693,0.0524]*1.0e-04;
        r_v = vertex_x1y1z1(index_list(i),:);%[0.0876,0.1668,0.0552]*1.0e-04;


        %equation for normal line
        syms t
        n = r0 - t*subs(fgrad,r,r0).'; % normal line
        mag_ = r_v-r0;
        mag_ = sqrt((mag_(1)^2)+(mag_(2)^2)+(mag_(3)^2));
        

        fplot3(n(1),n(2),n(3),[0 0.05],'b-','LineWidth',4);

        syms n(t)
        n(t) = r0 - t*subs(fgrad,r,r0).'; % normal line

        V = [r0;r_v];
        new = [n(0); n(1)];
        new = new/norm(new); %need to normalize
        
        %rotate about normal line
        
        final_temp = V-2*dot(V,new).*new;
        final = [n(0);final_temp(2,:)];
        % final=final_temp;

%         A = r0;
%         B = r_v; %starting line
%         C = final(2,:); %normal line
% 
%         S1 = B - A;
%         S2 = C - A;
%         Theta = atan2(norm(cross(S1, S2)), dot(S1, S2));

%         test = S2*sind(Theta);
        %find distance between point and normal line
        plot3([r0(1); final(2,1)],[r0(2); final(2,2)],[r0(3); final(2,3)],'g', 'LineWidth',3);
%         plot3([r0(1); test(1)],[r0(2); test(2)],[r0(3); test(3)],'g', 'LineWidth',3)


%         A = vertex_x0y0z0(index_list(i),:);
%         B = vertex_x1y1z1(index_list(i),:);
%         x = [V(1,1);V(2,1)];
%         y = [V(1,2);V(2,2)];
%         z = [V(1,3);V(2,3)];

%         plot3([V(1,1),P(2,1)],[V(1,2),P(2,2)],[V(1,),P(2,3)],'b.-','LineWidth',1)

%         legend("Original","Plane");%,"Normal")%,"Rotated")

        x1(index_list(i)) = final(2,1);%normal(1);
        y1(index_list(i)) = final(2,2);%normal(2);
        z1(index_list(i)) = final(2,3);%normal(3);
%         
%         
        
        result = [x1;y1;z1];
    end

end

% 
% %TODO:
% %Should probably vectorise this
% 
% function result = mirror_trajectory(previous_xyz,current_dxdydz, index_list)
%     
%     x0 = previous_xyz(1,:);
%     y0 = previous_xyz(2,:);
%     z0 = previous_xyz(3,:);
%     
%     x1 = x0 + current_dxdydz(1,:);
%     y1 = y0 + current_dxdydz(2,:);
%     z1 = z0 + current_dxdydz(3,:);
%     
%     result = [x1;y1;z1];
%     
%     vertex_x0y0z0 = [x0',y0',z0']; 
%     vertex_x1y1z1 = [x1',y1',z1'];
%     
% %     figure;
%     for i =1:length(index_list)
%         
%         V = [vertex_x0y0z0(index_list(i),:);...
%              vertex_x1y1z1(index_list(i),:)];
%          
% %         plot3(V(:,1),V(:,2),V(:,3),'k.-','MarkerSize',25,'Color', rand(1,3)); 
% %         xlabel('x'), ylabel('y'),zlabel('z');
% %         hold on
%         
%         A = vertex_x0y0z0(index_list(i),:);
%         B = vertex_x1y1z1(index_list(i),:);
%         x = [V(1,1);V(2,1)];
%         y = [V(1,2);V(2,2)];
%         z = [V(1,3);V(2,3)];
%         
%         normal = ([mean(x),mean(y),mean(z)] + null(A-B)');
%         normal = (normal(:,:)./sqrt(sum(normal(:,:).*normal(:,:))));%*1e-;%
%         
%         
% %         plot3([B(1),normal(1)],[B(2),normal(2)],[B(3),normal(3)],'r.-','LineWidth',1)
% %         legend("Original","Rotated")
% 
%         x1(index_list(i)) = normal(1);
%         y1(index_list(i)) = normal(2);
%         z1(index_list(i)) = normal(3);
%         result = [x1;y1;z1];
%     end
% 
% end
