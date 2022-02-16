%% Simulate Random Walk for n Particle that Starts Randomly in 3D Sphere
clear all; clc; close all; % clean up

num_particle      = 1000; 
START_TIME        = 0; %sec
STOP_TIME         = 0.05; %sec
movements_per_sec = 1000000;
dt = (STOP_TIME-START_TIME)/movements_per_sec;
t = movements_per_sec*dt;
numberOfSteps = (STOP_TIME-START_TIME)*movements_per_sec;

radii = (7*1e-6);
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
    sd__ = 0.58;%sqrt(STOP_TIME/numberOfSteps);

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
%                 XYZ = CheckCircleBounds(radii(radius),test_x,test_y,test_z,rand_x_steps(step,idx),rand_y_steps(step,idx),rand_z_steps(step,idx));
%                 xCoords(step, idx) =  XYZ(1,:);
%                 yCoords(step, idx) =  XYZ(2,:);
%                 zCoords(step, idx) =  XYZ(3,:);
                xCoords(step, idx) =  test_x;
                yCoords(step, idx) =  test_y;
                zCoords(step, idx) =  test_z;

                
            end %end of idx
        end


    %% Plotting avg walk for particles

    % figure; hold on
    
%     for particle = 1:num_particle
%         plot3(xCoords(:,particle),yCoords(:,particle),zCoords(:,particle),'Color', rand(1,3), 'MarkerSize', 9);
%         hold on; 
%     end
%     xlabel('x'), ylabel('y'),zlabel('z');
%     hold off
%     pause
    
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
    D_vec(1,radius) = MSD(radius)/(2*n*t)*(1/(STOP_TIME-START_TIME));


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