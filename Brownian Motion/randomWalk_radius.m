%% Simulate Random Walk for n Particle that Starts at Origin in 2D
clear all; clc; close all; % clean up

num_particle      = 20; 
START_TIME        = 0; %sec
STOP_TIME         = 200; %sec
movements_per_sec = 5;
numberOfSteps     = (STOP_TIME-START_TIME)*movements_per_sec;

radius = 1.5;


% vectors for random steps
rand_x_steps = (rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-0.5);
rand_y_steps = (rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-0.5);

extra_rand_x_steps = (rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-0.5);
extra_rand_y_steps = (rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-0.5);


xCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle); %particle start loc is assume 0,0
yCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle);

% DISK
x_center = 0; y_center = 0;

% loop through all steps
% TO DO:
    % - set a radius
for step = 2 : numberOfSteps
    for idx = 1:num_particle
            % x-loc
            
        test_x = xCoords(step-1, idx) + rand_x_steps(step,idx);
        test_y = yCoords(step-1, idx) + rand_y_steps(step,idx);

        % This section checks if any particles exist at the same place at
        % the same time.
        if idx > 1
            %checker_vec will be set to 1 where they match
            for checker = 1:(idx-1) %this section will probably need extention
                % - two particles cannot exist at same point  
                if ( (test_x == xCoords(step-1, checker) ) && (test_y == yCoords(step-1, checker)) )
                    disp("MATCH");
                    % x-loc
                    test_x = xCoords(step-1, idx) + extra_rand_x_steps(step,idx);%x1
                    % y-loc
                    test_y = yCoords(step-1, idx) + extra_rand_y_steps(step,idx);%y1
                end
            end %end of checker
        end %end of idx > 1

        % This section checks the set radius bounds
        CHECKER = CheckCircleBounds(x_center,test_x,y_center,test_y,radius);
        if (CHECKER == 1 || CHECKER == 2)
            % x-loc
            xCoords(step, idx) = xCoords(step-1, idx) + rand_x_steps(step,idx);%x1
            % y-loc
            yCoords(step, idx) = yCoords(step-1, idx) + rand_y_steps(step,idx);%y1

        elseif (CHECKER == 3)
            % is outside the circle
            % make the value opposite
            xCoords(step, idx) = xCoords(step-1, idx) - rand_x_steps(step,idx);%x1
            yCoords(step, idx) = yCoords(step-1, idx) - rand_y_steps(step,idx);%y1
        end            
    end %end of idx
end


%% Plotting walk for two particles

figure; hold on

for particle = 1:num_particle
    plot(xCoords(:,particle),yCoords(:,particle),'Color', rand(1,3), 'MarkerSize', 9);
    hold on; 
end
line(xlim, [0,0], 'Color', 'k', 'LineWidth', 1);
line([0,0], ylim, 'Color', 'k', 'LineWidth', 1);



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
 
function result = CheckCircleBounds(x0,x,y0,y,r)

    if ( sqrt((abs(x)-x0)^2 + (abs(y)-y0)^2) < r)
        result = 1; % (x,y) is inside the circle
    elseif ( sqrt((abs(x)-x0)^2 + (abs(y)-y0)^2) == r)
        result = 2; % (x,y) is on the circle
    elseif ( sqrt((abs(x)-x0)^2 + (abs(y)-y0)^2) > r)
        result = 3; % (x,y) is outside the circle
    end
    
end