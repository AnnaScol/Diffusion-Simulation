%% Simulate Random Walk for 1 Particle that Starts at Origin in 2D
clear all, clc


START_TIME        = 0; %sec
STOP_TIME         = 10; %sec
movements_per_sec = 5;
numberOfSteps = (STOP_TIME-START_TIME)*movements_per_sec;

rand_x_steps = 0.02*(rand((STOP_TIME-START_TIME)*movements_per_sec,1)-1);
rand_y_steps = 0.02*(rand((STOP_TIME-START_TIME)*movements_per_sec,1)-1);

particle_location = zeros(numberOfSteps,2);

xCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,1); %particle start loc is assume 0,0
yCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,1);


% loop through all steps
for step = 2 : numberOfSteps
    
	particle_location(step, 1) = particle_location(step-1, 1) + rand_x_steps(step);	
	particle_location(step, 2) = particle_location(step-1, 2) + rand_y_steps(step);
	% Now plot the walk so far.
	xCoords(step) = particle_location(step, 1);
	yCoords(step) = particle_location(step, 2);

end

%% Plot walked path
figure;
plot(xCoords(1),yCoords(1),'go-','LineWidth',2); 
hold on; 
plot(xCoords,yCoords,'bo-');
plot(xCoords(end),yCoords(end),'ro-','LineWidth',2);
xlim([min(xCoords)-0.005 max(xCoords)+0.005]), ylim([min(yCoords)-0.005 max(yCoords)+0.005]);
legend('Start Point (0,0)','Walked Path','End Point')
