%% Simulate Random Walk for n Particle that Starts at Origin in 3D
clear all; clc; close all; % clean up

num_particle      = 20; 
START_TIME        = 0; %sec
STOP_TIME         = 500; %sec
movements_per_sec = 8;
numberOfSteps     = (STOP_TIME-START_TIME)*movements_per_sec;

radius = 10;


% vectors for random steps
rand_x_steps = (rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-0.5);
rand_y_steps = (rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-0.5);
rand_z_steps = (rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-0.5);

% vectors for random steps if the new steps match
extra_rand_x_steps = (rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-0.5);
extra_rand_y_steps = (rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-0.5);
extra_rand_z_steps = (rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-0.5);

% vectors for finished random step placement for plotting & analysis
xCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle); %particle start loc is assume 0,0,0
yCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle);
zCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle);


% Sphere center
x_center = 0; y_center = 0; z_center = 0;

% Generate Sphere with hole
N=100;  MainSphereOrigin = [x_center,y_center,z_center];
CutOutRadius = 5; 
MainSphereRadius = 10;
[sphere_X,sphere_Y,sphere_Z,cutout_disk] = SphereHoles(MainSphereRadius,CutOutRadius,MainSphereOrigin,N);
cutout_idx = find(cutout_disk==1);
sphere_Z(cutout_idx) = 0;
%%
% loop through all steps
% TO DO:
    % - set a radius
for step = 2 : numberOfSteps
    
    for idx = 1:num_particle
        test_x = xCoords(step-1, idx) + rand_x_steps(step,idx);
        test_y = yCoords(step-1, idx) + rand_y_steps(step,idx);
        test_z = zCoords(step-1, idx) + rand_z_steps(step,idx);
        
        % This section checks the set radius bounds
        CHECKER = CheckSphereBounds([x_center,y_center,z_center],test_x,test_y,test_z,radius);
        
        if (CheckCutOut(MainSphereOrigin,test_x,test_y,test_z, sphere_X,sphere_Y,sphere_Z,cutout_idx) && CHECKER==3)
            %let it keep increasing as it is pasing the hole
            CHECKER = 1;
        end

        if (CHECKER == 1 || CHECKER == 2)
            xCoords(step, idx) = xCoords(step-1, idx) + rand_x_steps(step,idx);%x
            yCoords(step, idx) = yCoords(step-1, idx) + rand_y_steps(step,idx);%y
            zCoords(step, idx) = zCoords(step-1, idx) + rand_z_steps(step,idx);%z

        elseif (CHECKER == 3)
            % is outside the circle
            % make the value opposite
            xCoords(step, idx) = xCoords(step-1, idx) - rand_x_steps(step,idx);%x
            yCoords(step, idx) = yCoords(step-1, idx) - rand_y_steps(step,idx);%y
            zCoords(step, idx) = zCoords(step-1, idx) - rand_z_steps(step,idx);%z
        end            
    end %end of idx
    
end


%% Plotting walk for two particles

figure; 

for particle = 1:num_particle
    plot3(xCoords(:,particle),yCoords(:,particle),zCoords(:,particle));
    hold on; 
end

%plot line axis
line_ = -15:0.01:15;
other_coords = zeros(1,length(line_));
plot3(line_,other_coords,other_coords,'Color', 'k', 'LineWidth', 1);
plot3(other_coords,line_,other_coords,'Color', 'k', 'LineWidth', 1);
plot3(other_coords,other_coords,line_,'Color', 'k', 'LineWidth', 1);

mesh(sphere_X,sphere_Y,sphere_Z,'edgealpha',0.3,'facealpha',0.5)

hold off

%% Finding average distances from origin 
distancesFromOrigin = hypot(xCoords(end,:), yCoords(end,:));
figure;
histObject = histogram(distancesFromOrigin, 25);
grid on;
caption = sprintf('Distribution of %d Final Distances', num_particle);
title(caption); xlabel('Distance'); ylabel('Count');

averageFinalX = mean(xCoords(end,:));
averageFinalY = mean(yCoords(end,:));



%%
%check if value is within circle
% 1) if (x-x0)^2 + (y-y0)^2 < r^2, the point (x,y) is inside the circle,
% 2) if (x-x0)^2 + (y-y0)^2 == r^2, the point (x,y) is on the circle, and
% 3) if (x-x0)^2 + (y-y0)^2 > r^2, the point (x,y) is outside the circle.
 
function result = CheckSphereBounds(origin,x,y,z,r)

    mag = sqrt((abs(x)-origin(1))^2 + (abs(y)-origin(2))^2 + (abs(z)-origin(3))^2);
    if ( mag < r)
        result = 1; % (x,y) is inside the circle
    elseif ( mag == r)
        result = 2; % (x,y) is on the circle
    elseif ( mag > r)
        result = 3; % (x,y) is outside the circle
    end
    
end


function result = CheckCutOut(origin,x,y,z, sphere_X,sphere_Y,sphere_Z,cut_out_idx)
%cutoutind is [x,y]
%just assume start at origin

    d = (x-sphere_X).^2+(y-sphere_Y).^2+(z-sphere_Z).^2; %// compute squared distances
%     if max(d,[],'omitnan') > 10^2
%         disp("beep")
%     end
    result = d;   
    cutout_test_case = x^2+y^2+z^2;

    %if the magntidue is equal to the radius then we can check if its
    %passing point is the hole, where the hole indicies are at cut_out_idx
    %defined as a vector
%     if ((cutout_test_case >= 110)) 
%         if (d(cut_out_idx) <= 3)
%             result = 1;
%             disp("beep");
%         end 
%     else
%         result = 0;
%     end

    if ((y >= 10) && (abs(x) < 3))
%         if (d(cut_out_idx) <= 3)
            result = 1;
%             disp("beep");
%         end 
    else
        result = 0;
    end
    
    
    
         
end