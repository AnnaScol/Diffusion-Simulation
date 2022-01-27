%% Simulate Random Walk for n Particle that Starts at Origin in 3D
clear all; clc; close all; % clean up

num_particle      = 40; 
START_TIME        = 0; %sec
STOP_TIME         = 700; %sec
movements_per_sec = 4;
numberOfSteps     = (STOP_TIME-START_TIME)*movements_per_sec;

radius = 8;
step_count = (STOP_TIME-START_TIME)*movements_per_sec;

% vectors for random steps
rand_x_steps = (rand(step_count,num_particle)-0.5);
rand_y_steps = (rand(step_count,num_particle)-0.5);
rand_z_steps = (rand(step_count,num_particle)-0.5);

% vectors for random steps if the new steps match
extra_rand_x_steps = (rand(step_count,num_particle)-0.5);
extra_rand_y_steps = (rand(step_count,num_particle)-0.5);
extra_rand_z_steps = (rand(step_count,num_particle)-0.5);

% vectors for finished random step placement for plotting & analysis
xCoords = zeros(step_count,num_particle); %particle start loc is assume 0,0,0
yCoords = zeros(step_count,num_particle);
zCoords = zeros(step_count,num_particle);


% Sphere center
x_center = 0; y_center = 0; z_center = 0;

% Generate Sphere with hole
N=80;  MainSphereOrigin = [x_center,y_center,z_center];
MainSphereRadius = 8;
CuttOutRadius = [4 4 4];
CuttOutCenter = [10 30 -30];

[sphere_X,sphere_Y,sphere_Z,cutout_disk,cutout_idx] = CutOutSphere(MainSphereRadius,CuttOutRadius,CuttOutCenter,N);

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
        
        if (CHECKER==3 && CheckCutOut(MainSphereOrigin,test_x,test_y,test_z, sphere_X,sphere_Y,sphere_Z,cutout_disk,cutout_idx) )
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

for i = 1:length(CuttOutRadius)
    sphere_X(cutout_idx(i,(find(cutout_idx(i,:)>1)))) = nan;
    sphere_Y(cutout_idx(i,(find(cutout_idx(i,:)>1)))) = nan;
    sphere_Z(cutout_idx(i,(find(cutout_idx(i,:)>1)))) = nan; 
end

mesh(sphere_X,sphere_Y,sphere_Z,'edgealpha',0.8,'facealpha',0.8)
xlabel('x')
ylabel('y')
zlabel('z')
hold off

%% Finding average distances from origin 
distancesFromOrigin = hypot(xCoords(end,:), yCoords(end,:));
figure;
histObject = histogram(distancesFromOrigin, 40);
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
%need to let particle out if it is within the cut out region
function result = CheckCutOut(origin,x,y,z, sphere_X,sphere_Y,sphere_Z,cutout_disk,cutout_idx)
%just assume start at origin
    result = 0;

    %try testing each of the x,y,x and x,y,z sphere values in a loop
    for i = 1:size(cutout_idx,1)
        
% cutout_idx(i,(find(cutout_idx(i,:)>1)))
        
        S_index = cutout_idx(i,find(cutout_idx(i,:)~=0));
        new_S_index = 1:length(find(squeeze(cutout_disk(i,:,:))==1));
        
        sphere_X_val(i,new_S_index) = sphere_X(S_index);
        sphere_Y_val(i,new_S_index) = sphere_Y(S_index);
        sphere_Z_val(i,new_S_index) = sphere_Z(S_index);
        
    end
    
    
    for num_cutouts = 1:size(cutout_disk,1)
        
        for i = 1:length(find(sphere_Z_val(num_cutouts,:)~=0))
            
            Sx = sphere_X_val(num_cutouts,i);
            Sy = sphere_Y_val(num_cutouts,i);
            Sz = sphere_Z_val(num_cutouts,i);
            
            %check if it is meant to be negative or positive in the axis limits
            if ( (Sx >= 0) && (Sy >= 0) && (Sz >= 0) )
                
                if ( (Sx <= x) && (Sy <= y) && (Sz <= z) )
%                     fprintf("Logic %d, Sx = %.2f, Sy = %.2f, Sz = %.2f, x=%.2f y=%.2f z=%.2f\n",((Sx <= x) && (Sy <= y) && (Sz <= z) ),Sx,Sy,Sz,x,y,z)
                    result = 1;
                end
                
            elseif ( (Sx >= 0) && (Sy >= 0) && (Sz < 0) )
                
                if ( (Sx <= x) && (Sy <= y) && (Sz >= z) )
                    result = 1;
                end

            elseif ( (Sx >= 0) && (Sy < 0) && (Sz >= 0) )
                
                if ( (Sx <= x) && (Sy >= y) && (Sz <= z) )
                    result = 1;
                end
                
            elseif ( (Sx < 0) && (Sy >= 0) && (Sz >= 0) )
                
                if ( (Sx >= x) && (Sy <= y) && (Sz <= z) )
                    result = 1;
                end      
               
            elseif ( (Sx >= 0) && (Sy < 0) && (Sz < 0) )
                
                if ( (Sx <= x) && (Sy >= y) && (Sz >= z) )
                    result = 1;
                end   
                
            elseif ( (Sx < 0) && (Sy < 0) && (Sz >= 0) )
                
                if ( (Sx  >= x) && (Sy >= y) && (Sz <= z) )
                    result = 1;
                end 
                
            elseif ( (Sx < 0) && (Sy >= 0) && (Sz < 0) )
                
                if ( (Sx  >= x) && (Sy <= y) && (Sz >= z) )
                    result = 1;
                end 
            
            elseif ( (Sx < 0) && (Sy < 0) && (Sz < 0) )
                
                if ( (Sx >= x) && (Sy >= y) && (Sz >= z) )
                    result = 1;
                end  
                
            end %end of if,elseif,..... statement
        end %end of all index checks
    end %end of num_counts 

end


