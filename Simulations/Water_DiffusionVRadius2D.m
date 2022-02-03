%% Simulate Random Walk for n Particle that Starts at Origin in 2D
clear all; clc; close all; % clean up

num_particle      = 1000; 
START_TIME        = 0; %sec
STOP_TIME         = 1; %sec
movements_per_sec = 10000;
dt = (STOP_TIME-START_TIME)/movements_per_sec;
t = movements_per_sec*dt;
numberOfSteps = (STOP_TIME-START_TIME)*movements_per_sec;

% radii = ((50)*1e-6);
radii = (linspace(0.001,1,100)*1e-6);

MSD = zeros(1,length(radii));
D_vec = zeros(1,length(radii));

D = 3e-9; %free water
n = 2;

xCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle); %particle start loc is assume 0,0
yCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle); %particle start loc is assume 0,0

time = (1:size(xCoords,1))*dt;

nTests = 100;
resultStorage_x = zeros(nTests,(STOP_TIME-START_TIME)*movements_per_sec,num_particle);
resultStorage_y = zeros(nTests,(STOP_TIME-START_TIME)*movements_per_sec,num_particle);
%0.0320
for radius = 1:length(radii)
    fprintf("RADIUS %d/%d\n",radius,length(radii));

    
        rand_x_steps = ((0.0059/1).*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle))*sqrt(2*n*D*dt);
        rand_y_steps = ((0.0059/1).*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle))*sqrt(2*n*D*dt);
        
        test_x = rand_x_steps(1,:);
        test_y = rand_y_steps(1,:);
        
%         XY = CheckCircleBounds(radii(radius),test_x,test_y,rand_x_steps(1,idx),rand_y_steps(1,idx));
%         xCoords(step, idx) =  XY(1,:);
%         yCoords(step, idx) =  XY(2,:);
        
%         rand_x_steps = ((0.0323/2).*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle))*sqrt(2*n*D*dt);
%         rand_y_steps = ((0.0323/2).*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle))*sqrt(2*n*D*dt);
%         
                
        for step = 2:numberOfSteps
            for idx = 1:num_particle

                test_x = xCoords(step-1, idx) + rand_x_steps(step,idx);
                test_y = yCoords(step-1, idx) + rand_y_steps(step,idx);

                % - two particles cannot exist at same point    
                %Radius Bounds       
                XY = CheckCircleBounds(radii(radius),test_x,test_y,rand_x_steps(step,idx),rand_y_steps(step,idx));
                xCoords(step, idx) =  XY(1,:);
                yCoords(step, idx) =  XY(2,:);
                
%                 if (CHECKER == 1)
%                     % x-loc
%                     xCoords(step, idx) = xCoords(step-1, idx) + rand_x_steps(step,idx);%x1
%                     % y-loc
%                     yCoords(step, idx) = yCoords(step-1, idx) + rand_y_steps(step,idx);%y1
% 
%                 elseif (CHECKER == 3)
%                     % is outside the circle
%                     % make the value opposite
%                     xCoords(step, idx) = xCoords(step-1, idx) - rand_x_steps(step,idx);%x1
%                     yCoords(step, idx) = yCoords(step-1, idx) - rand_y_steps(step,idx);%y1
%                 end   
                
            end %end of idx
        end


    %% Plotting avg walk for particles

    % figure; hold on
    
    for particle = 1:num_particle
        plot(xCoords(:,particle),yCoords(:,particle),'Color', rand(1,3), 'MarkerSize', 9);
        hold on; 
    end
    line(xlim, [0,0], 'Color', 'k', 'LineWidth', 1);
    line([0,0], ylim, 'Color', 'k', 'LineWidth', 1);
    hold off
%     pause
    
% 
    % Calculate the distance from the origin.
    distancesFromOrigin_x = xCoords(end,:);
    distancesFromOrigin_y = yCoords(end,:);
    magnitudeFromOrigin = hypot(distancesFromOrigin_x,distancesFromOrigin_y);
    
% % 
%     figure;
%     histObject = histogram(magnitudeFromOrigin, 25);
%     grid on;
%     caption = sprintf('Distribution of %d Final Distances', num_particle);
%     title(caption);
%     xlabel('Distance');ylabel('Count');
%     
    MSD(1,radius) = mean(abs(magnitudeFromOrigin).^2)
    D_vec(1,radius) = (MSD(radius)/(2*n*dt))*(1/t)
end
%%
figure
plot(radii,D_vec,'.-');
hold on;
xlabel('Radius Size (m)');ylabel('D(m^{2}/s)')
line([0,radii(end)], [D, D], 'Color', 'k', 'LineWidth', 1);
title("Diffusivity varying over Radii Constraints");
legend('D-Measurements', 'Free water D = 3e-9');
saveas(gcf,'1000Spin_10000ts_100p.fig')

%%
%check if value is within circle
% 1) if (x-x0)^2 + (y-y0)^2 < r^2, the point (x,y) is inside the circle,
% 2) if (x-x0)^2 + (y-y0)^2 == r^2, the point (x,y) is on the circle, and
% 3) if (x-x0)^2 + (y-y0)^2 > r^2, the point (x,y) is outside the circle.
%always at origin
function result = CheckCircleBounds(r,x,y,dx,dy)
    test = (sqrt((abs(x))^2 + (abs(y))^2) < r);
    if ( test == 1)% (x,y) is inside the circle
        result = [x;y]; 
    elseif ( test == 0) % (x,y) is outside or equal to the circle
        result = mirror_trajectory([x,y],[dx dy]);
    end
    
end

function result = mirror_trajectory(current_xy,current_dxdy)
    
    x0 = current_xy(1)-current_dxdy(1);
    y0 = current_xy(2)-current_dxdy(2);
    
    x1 = current_dxdy(1);
    y1 = current_dxdy(2);
    
    result = [x1;y1];
    
    vertex_x0y0 = [x0,y0]; 
    vertex_x1y1 = [x1,y1];
    
    V = [vertex_x0y0;...
         vertex_x1y1];
     
%     figure; 
%     plot(V(:,1),V(:,2),'k.-','MarkerSize',25,'Color', rand(1,3)); 
%     xlabel('x'), ylabel('y')
%     hold on

    A = vertex_x0y0;
    B = vertex_x1y1;
    x = [V(1,1);V(2,1)];
    y = [V(1,2);V(2,2)];
    
    normal = ([mean(x),mean(y)] + null(A-B)');
    normal = [B(1), B(2);...
                normal];

    r = V-2*(dot(V,normal))*normal;
%     normal = ([mean(x),mean(y)] + null(A-B)');
%     normal = (normal(:,:)./sqrt(sum(normal(:,:).*normal(:,:))))*1e-6;%

% 
%     plot([B(1),normal(1)*1e-6],[B(2),normal(2)*1e-6],'r.-','LineWidth',1)
%     plot([B(1),r(2)],[B(2),r(2)],'b.-','LineWidth',1)
% 
%     legend("Original","Rotated")

    result = [r(1);r(2)];
end