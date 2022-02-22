%% Simulate Random Walk for 2 Particle that Starts at Origin in 2D
clear all; clc; close all; % clean up

num_particle      = 2; 
START_TIME        = 0; %sec
STOP_TIME         = 10; %sec
movements_per_sec = 1;
numberOfSteps = (STOP_TIME-START_TIME)*movements_per_sec;

% vectors for random steps
rand_x_steps = 0.02*(rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-1);
rand_y_steps = 0.02*(rand((STOP_TIME-START_TIME)*movements_per_sec,num_particle)-1);

xCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle); %particle start loc is assume 0,0
yCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle);


% loop through all steps
for step = 2 : numberOfSteps
    for idx = 1:num_particle
        % x-loc
        xCoords(step, idx) = xCoords(step-1, idx) + rand_x_steps(step,idx);%x1
        
        % y-loc
        yCoords(step, idx) = yCoords(step-1, idx) + rand_y_steps(step,idx);%y1
    end 
end


%% Plotting walk for two particles
colour = ['b','k'];

figure; hold on

for particle = 1:num_particle
    plot(xCoords(1,particle),xCoords(1,particle),'g*-','LineWidth',2); 
    plot(xCoords(:,particle),yCoords(:,particle),sprintf('%co-',colour(particle)));
    plot(xCoords(end,particle),yCoords(end,particle),'ro-','LineWidth',2);
    hold on; 
end
hold off

