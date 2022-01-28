%% Simulate Random Walk for n Particle that Starts at Origin in 2D
clear all; clc; close all; % clean up

num_particle      = 20; 
START_TIME        = 0; %sec
STOP_TIME         = 5; %sec
movements_per_sec = 1;
numberOfSteps = (STOP_TIME-START_TIME)*movements_per_sec;

% vectors for random steps
rand_x_steps = (randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle));
rand_y_steps = (randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle));

xCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec+20,num_particle); %particle start loc is assume 0,0
yCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec+20,num_particle);


% loop through all steps
% TO DO:
    % - set a radius
for step = 2 : numberOfSteps
    for idx = 1:num_particle
            % x-loc
            xCoords(step, idx) = xCoords(step-1, idx) + rand_x_steps(step,idx);%x1
            % y-loc
            yCoords(step, idx) = yCoords(step-1, idx) + rand_y_steps(step,idx);%y1
        
        % - two particles cannot exist at same point    
        if idx > 1
            
            %checker_vec will be set to 1 where they match
            checker_vec = zeros(1,idx);
            all_zeros =  zeros(1,idx);
            for checker = 1:idx
                if (xCoords(step, idx) == xCoords(step-1, idx-1) && yCoords(step, idx) == yCoords(step-1, idx-1))
                % x-loc
                xCoords(step, idx) = xCoords(step-1, idx) + rand_x_steps(step+1,idx);%x1
                % y-loc
                yCoords(step, idx) = yCoords(step-1, idx) + rand_y_steps(step+1,idx);%y1
                end
                
            end %end of checker
        end %end of idx > 1
    end %end of idx
end


%% Plotting walk for two particles
colour = ['b','k'];

figure; hold on

for particle = 1:num_particle
    plot(xCoords(:,particle),yCoords(:,particle),'Color', rand(1,3), 'MarkerSize', 9);
    hold on; 
end
line(xlim, [0,0], 'Color', 'k', 'LineWidth', 1);
line([0,0], ylim, 'Color', 'k', 'LineWidth', 1);
hold off
