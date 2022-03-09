%% Simulate Random Walk for n Particle that Starts at Origin in 2D
clear all; clc; close all; % clean up

num_particle      = 1000; 
START_TIME        = 0; %sec
STOP_TIME         = 0.06; %sec ldelta+sdelta+sdelta with 1/t for final d std is 0.0072/3.18, without it is
movements_per_sec = 100000;
dt = (STOP_TIME-START_TIME)/movements_per_sec;
t = movements_per_sec*dt;
numberOfSteps = (STOP_TIME-START_TIME)*movements_per_sec;

% radii = ((20)*1e-6);
% radii = (linspace(1,8,20)*1e-6);
% radii = (linspace(0.01,10,100)*1e-6);
radii = ((10)*1e-6);


MSD = zeros(1,length(radii));
D_vec = zeros(1,length(radii));

D = 3e-9; %free water
n = 2;

xCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle); %particle start loc is assume 0,0
yCoords = zeros((STOP_TIME-START_TIME)*movements_per_sec,num_particle); %particle start loc is assume 0,0

time = (1:size(xCoords,1))*dt;

%0.0320
for radius = 1:length(radii)
    fprintf("RADIUS %d/%d\n",radius,length(radii));

    
        rand_x_steps = ((0.0095).*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle))*sqrt(2*n*D*dt);
        rand_y_steps = ((0.0095).*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle))*sqrt(2*n*D*dt);
        
%         r = (radii(radius)-1.0e-7) * sqrt(rand(1,num_particle));
%         theta = rand(1,num_particle)*2*pi;
%         x = zeros(1,num_particle) + r.*cos(theta);
%         y = zeros(1,num_particle) + r.*sin(theta);
%         
%         
%         xCoords(1,:) = x;
%         yCoords(1,:) = y;
        xCoords(1,:) = zeros(1,num_particle);
        yCoords(1,:) = zeros(1,num_particle);
        plot(xCoords(1,:),yCoords(1,:),'.')
        
%         for p = 1:num_particle
%             XY = CheckCircleBounds(radii(radius),test_x(1,p),test_y(1,p),rand_x_steps(1,p),rand_y_steps(1,p));
%             xCoords(1, p) =  XY(1,:);
%             yCoords(1, p) =  XY(2,:);
%         end
%         
%         rand_x_steps = ((0.0323/2).*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle))*sqrt(2*n*D*dt);
%         rand_y_steps = ((0.0323/2).*randn((STOP_TIME-START_TIME)*movements_per_sec,num_particle))*sqrt(2*n*D*dt);
%         
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

    magnitudeFromOrigin = sqrt((abs(distancesFromOrigin_x - xCoords(1,:)).^2) ...
                             + (abs(distancesFromOrigin_y - yCoords(1,:)).^2));

    
%     magnitudeFromOrigin = hypot(distancesFromOrigin_x,distancesFromOrigin_y);
    
% % 
%     figure;
%     histObject = histogram(magnitudeFromOrigin, 25);
%     grid on;
%     caption = sprintf('Distribution of %d Final Distances', num_particle);
%     title(caption);
%     xlabel('Distance');ylabel('Count');
%     
    MSD(1,radius) = mean(abs(magnitudeFromOrigin).^2)
    D_vec(1,radius) = (MSD(radius)/(2*n*dt))%*(1/t)
end
%%
figure
plot(radii,D_vec,'.-');
hold on;
xlabel('Radius Size (m)');ylabel('D(m^{2}/s)')
% line([0,radii(end)], [D, D], 'Color', 'k', 'LineWidth', 1);
title("Diffusivity varying over Radii Constraints");
legend('D-Measurements', 'Free water D = 3e-9');
saveas(gcf,'uniform_start.fig')

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
%         result = mirror_trajectory([x,y],[dx dy]);
        result = [x-dx;y-dy];
    end
    
end

function result = mirror_trajectory(current_xy,current_dxdy)
    
    x0 = current_xy(1)-current_dxdy(1);
    y0 = current_xy(2)-current_dxdy(2);
    
    x1 = current_xy(1);
    y1 = current_xy(2);
    
    result = [x1;y1];
    
    vertex_x0y0 = [x0,y0]; 
    vertex_x1y1 = [x1,y1];
    
    V = [vertex_x0y0;...
         vertex_x1y1];

    A = vertex_x0y0;
    B = vertex_x1y1;
    x = [V(1,1);V(2,1)];
    y = [V(1,2);V(2,2)];
    
    normal = ([mean(x),mean(y)] + null(A-B)');
    normal = [B(1), B(2);...
                normal];
    orig_mag = (sqrt((abs(x1-x0))^2 + (abs(y1-y0))^2));
    
    normal = normal*orig_mag;
    r = V-2*(dot(V,normal))*normal;
    r = r*1e-20;

    
%     theta =  10;

% %     theta =  -atan2d(V(1,2)-normal(2,2),V(1,1)-normal(2,1))*1e-5;
% %     R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
% %     vR = V(1,:)*R;

%     figure; 
%     plot([x0 x1],[y0 y1],'k.-','MarkerSize',25,'Color', rand(1,3)); 
%     xlabel('x'), ylabel('y')
%     hold on
%     plot([normal(2,1) x1],[normal(2,2) y1],'r.-','LineWidth',1)
%     plot([(r(1,1)) x1],[(r(1,2)) y1],'b.-','LineWidth',1)
%     legend("Original","normal","Rot1")
% %     

    result = [vR(1);vR(2)];
end