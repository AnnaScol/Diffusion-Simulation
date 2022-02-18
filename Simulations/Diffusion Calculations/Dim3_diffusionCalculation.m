%% Simulate Random Walk for n Particle that Starts Randomly in 3D Sphere
clear all; clc; close all; % clean up

num_particle      = 500; 
START_TIME        = 0; %sec
STOP_TIME         = 0.05; %sec
movements_per_sec = 1000000;
dt = (STOP_TIME-START_TIME)/movements_per_sec;
t = movements_per_sec*dt;
numberOfSteps = (STOP_TIME-START_TIME)*movements_per_sec;

radii = (10*1e-6);
% radii = (linspace(1,15,30)*1e-6);

MSD = zeros(1,length(radii));
D_vec = zeros(1,length(radii));

D = 3e-9; %free water
n = 3;

xCoords = zeros(numberOfSteps,num_particle); %particle start loc is assume 0,0
yCoords = zeros(numberOfSteps,num_particle); %particle start loc is assume 0,0
zCoords = zeros(numberOfSteps,num_particle); %particle start loc is assume 0,0

for radius = 1:length(radii)
    fprintf("RADIUS %d/%d\n",radius,length(radii));
%     sd__ = 1.85;%sqrt(STOP_TIME/numberOfSteps);
    sd__ = 0.58;%sqrt(STOP_TIME/numberOfSteps); was 0.58

    rand_x_steps = sd__.*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle)*sqrt(2*n*D*dt);
    rand_y_steps = sd__.*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle)*sqrt(2*n*D*dt);
    rand_z_steps = sd__.*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle)*sqrt(2*n*D*dt);
            
%         point_in_sphere = StartingPos(radii(radius),num_particle);
    xCoords(1,:) = zeros(1,num_particle);
    yCoords(1,:) = zeros(1,num_particle);
    zCoords(1,:) = zeros(1,num_particle);
%             
%         point_in_sphere = StartingPos(radii(radius),num_particle);
%         xCoords(1,:) = point_in_sphere(1,:)';
%         yCoords(1,:) = point_in_sphere(2,:)';
%         zCoords(1,:) = point_in_sphere(3,:)';
        
        for step = 2:numberOfSteps
            for idx = 1:num_particle

                test_x = xCoords(step-1, idx) + rand_x_steps(step,idx);
                test_y = yCoords(step-1, idx) + rand_y_steps(step,idx);
                test_z = zCoords(step-1, idx) + rand_z_steps(step,idx);

                % - two particles cannot exist at same point    
                %Radius Bounds       
                XYZ = InBounds(radii(radius),test_x,test_y,test_z,rand_x_steps(step,idx),rand_y_steps(step,idx),rand_z_steps(step,idx));
                xCoords(step, idx) =  XYZ(1,:);
                yCoords(step, idx) =  XYZ(2,:);
                zCoords(step, idx) =  XYZ(3,:);
%                 xCoords(step, idx) =  test_x;
%                 yCoords(step, idx) =  test_y;
%                 zCoords(step, idx) =  test_z;

                
            end %end of idx
        end


    %% Plotting avg walk for particles

    figure; hold on
    
    for particle = 1:num_particle
        plot3(xCoords(:,particle),yCoords(:,particle),zCoords(:,particle),'Color', rand(1,3), 'MarkerSize', 9);
        hold on; 
    end
    xlabel('x'), ylabel('y'),zlabel('z');
    hold off
   
    
% 
    % Calculate the distance from the origin.
    distancesFromOrigin_x = xCoords(end,:)-xCoords(1,:);
    distancesFromOrigin_y = yCoords(end,:)-yCoords(1,:);
    distancesFromOrigin_z = zCoords(end,:)-zCoords(1,:);
    

    magnitudeFromOrigin = sqrt((abs(distancesFromOrigin_x)).^2 +...
                               (abs(distancesFromOrigin_y)).^2 +...
                               (abs(distancesFromOrigin_z)).^2);
                              
% % 
    figure;
    histObject = histogram(magnitudeFromOrigin, 25);
    grid on;
    caption = sprintf('Distribution of %d Final Distances', num_particle);
    title(caption);
    xlabel('Distance');ylabel('Count');
    
    
    dr_squared = (distancesFromOrigin_x).^2 +...
                 (distancesFromOrigin_y).^2 +...
                 (distancesFromOrigin_z).^2;
                      
                      
    MSD(1,radius) = mean(dr_squared);
%     D_vec(1,radius) = (MSD(radius)/(2*n*dt))*(1/t); WHAT IS CHANGED
%     D_vec(1,radius) = MSD(radius)/(2*n*t)*(1/(STOP_TIME-START_TIME));
    D_vec(1,radius) = MSD(radius)/(2*n*numberOfSteps*dt);%*(1/(STOP_TIME-START_TIME));



    disp(D_vec);
end
%%
figure
plot(radii,D_vec,'.-');
hold on;
xlabel('Radius Size (m)');ylabel('D(m^{2}/s)')
line([0,radii(end)], [D, D], 'Color', 'k', 'LineWidth', 1);
title("Diffusivity varying over Radii Constraints");
legend('D-Measurements', 'Free water D = 3e-9');

%%



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





%check if value is within circle
% 1) if (x-x0)^2 + (y-y0)^2 < r^2, the point (x,y) is inside the circle,
% 2) if (x-x0)^2 + (y-y0)^2 == r^2, the point (x,y) is on the circle, and
% 3) if (x-x0)^2 + (y-y0)^2 > r^2, the point (x,y) is outside the circle.
%always at origin
function result = CheckCircleBounds(r,x1,y1,z1,dx,dy,dz)

    result = [x1;y1;z1];
    test_axis_x = sqrt((abs(x1)).^2 + (abs(y1-dy)).^2 + (abs(z1-dz)).^2);
    test_axis_y = sqrt((abs(x1-dx)).^2 + (abs(y1)).^2 + (abs(z1-dy)).^2);
    test_axis_z = sqrt((abs(x1-dx)).^2 + ((abs(y1-dy))).^2 + (abs(z1)).^2);
    test_axis = [test_axis_x;test_axis_y;test_axis_z];
    test_res = (test_axis>=r);
    temp = max(test_res);
    temp = max(temp);
    
%     if (temp == 1)
%         
%         mirrored = mirror_trajectory([(x1-dx);(y1-dy);(z1-dz)], [dx;dy;dz]);
%         result = mirrored;
%     end

% 
%     test = (sqrt((abs(x))^2 + (abs(y))^2) < r);
%     if ( test == 1)% (x,y) is inside the circle
%         result = [x;y]; 
%     elseif ( test == 0) % (x,y) is outside or equal to the circle
%         result = mirror_trajectory([x,y,z],[dx dy dz]);
%     end
    
end



function result = mirror_trajectory(current_xyz,current_dxdydz)
    
    x0 = current_xyz(1)-current_dxdydz(1);
    y0 = current_xyz(2)-current_dxdydz(2);
    z0 = current_xyz(2)-current_dxdydz(2);
    
    x1 = current_dxdydz(1);
    y1 = current_dxdydz(2);
    z1 = current_dxdydz(2);
    
    result = [x1;y1;z1];
    
    vertex_x0y0z0 = [x0,y0,z0]; 
    vertex_x1y1z1 = [x1,y1,z1];
    
    V = [vertex_x0y0z0;...
         vertex_x1y1z1];

%     figure;
%     plot3(V(:,1),V(:,2),V(:,3),'k.-','MarkerSize',25,'Color', rand(1,3)); 
%     xlabel('x'), ylabel('y'),zlabel('z');
%     hold on

    A = vertex_x0y0z0;
    B = vertex_x1y1z1;
    x = [V(1,1);V(2,1)];
    y = [V(1,2);V(2,2)];
    z = [V(1,3);V(2,3)];
    
    normal = ([mean(x),mean(y),mean(z)] + null(A-B)');
    normal = [B(1), B(2), B(3);...
              normal(2,1),normal(2,2),normal(2,3)];

    r = V-2*(dot(V,normal)).*normal;

%     plot3([B(1),normal(1)*1e-6],[B(2),normal(2)*1e-6],[B(3),normal(3)*1e-6],'r.-','LineWidth',1)
%     plot3([B(1),r(2)],[B(2),r(2)],[B(3),r(3)],'b.-','LineWidth',1)
%     legend("Original","Normal","rotated")
% 
%     legend("Original","Rotated")

    result = [r(1);r(2);r(3)];
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