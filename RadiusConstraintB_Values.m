%% Spin Echo Line Scan
% TODO:
% Vectorize the mirror bounce 

clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 10*10^-6; 
gamma = 2*pi*42.577*10^6;


%Allocate the memory needed
nTimeSteps  = 7000*2*10^-5/dt;%70ms
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
adc         = zeros(1,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points

% %% Parameters %% %

% Diffusion Gradients
ldelta = 0.040; %ms
sdelta = 0.020; %ms
D = 1e-6;    %m^2/s

nSpins = 500;
G = ([0,5,7.5,10,12.5,15,17.5,20,22,25,30,35,40,45,50,55,60,65,70,75,80,100,125,150,175,200,300]*1e-3); %mT


nG = length(G);
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position

b = gamma^2*G.^2*sdelta^2*(ldelta-sdelta/3);

T1 = 1000*(10^-3);
T2 = 1000*(10^-3); 
TE = 65*(10^-3);

% Make size soma cell and brain cell sizes - vary 0.005-0.1mm
constraint_radii = ([1 5 10 15 20 25 35 1*1e6]*1e-6); 
% %%   %% %
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

% Generate and Display MRI Sequence
% %% RF %% %
% RF Excitation Waveform
pulsedurE = 0.002; % duration of the RF in s
rfStepsE = round(1:(pulsedurE/(dt)));
rfPulse(rfStepsE) = apodize_sinc_rf(length(rfStepsE),3,pi/2,dt); %B1+ in Tesla

% First diffusion gradient
diffusionGradient1_loc = round((1:(sdelta/dt))+ pulsedurE/dt);
gradAmp(3,diffusionGradient1_loc) =  G(1); %Z gradients in Tesla per meter

%RF Refocusing pules
pulsedurR = 0.002; % duration of the RF in s
rfStepsR = round(1:(pulsedurR/(dt)));
rfPulseR = apodize_sinc_rf(length(rfStepsR),3,pi,dt); %B1+ in Tesla
rfPulse(round(TE/2/dt) +length(rfStepsE)/2 + rfStepsR) = rfPulseR;

% diffusion pulse 2
diffusionGradient2_loc = round((1:(sdelta/dt)) + pulsedurE/dt + pulsedurR/dt + ldelta/dt + sdelta/dt);
gradAmp(3,diffusionGradient2_loc) =  G(1); %Z gradients in Tesla per meter

location = zeros(3,nTimeSteps);

%%% Plotting the data points %%%
Coords = zeros(length(constraint_radii),3,nSpins,nTimeSteps); %particle start loc is assume 0,0,0
% yCoords = zeros(length(constraint_radii),nSpins,nTimeSteps);
% zCoords = zeros(length(constraint_radii),nSpins,nTimeSteps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% PLOTTING %% %
% figure
% subplot(3,1,1); plot(time,rfPulse,'k-','LineWidth',2);title('RF Pulse'); 
% xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;
% 
% subplot(3,1,2); plot(time,phase(rfPulse),'k-','LineWidth',2);title('RF Pulse Phase'); 
% xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;
% 
% subplot(3,1,3); plot(time,gradAmp(3,:),'b-','LineWidth',2);title('Slice Select Gradient');
% xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;

% %% Get location matrix --> 3 x nSpins dimensions
SetsOfSpins=20;
store_final_vect = zeros(SetsOfSpins,length(b),length(constraint_radii));

%%
for set = 1:20 %20 sets of 500
    
    for radius_bounds = 1:length(constraint_radii)
            j = 1;
            fprintf("RADIUS %d",radius_bounds);
            %starting spin locations
            r = zeros(3,nSpins);
            Coords(radius_bounds,1,:,j) = r(1,:)';
            Coords(radius_bounds,2,:,j) = r(2,:)';
            Coords(radius_bounds,3,:,j) = r(3,:)';

            for j = 2:nTimeSteps   
                %%% Bounds checking %%%
                rnd = (-1 + 2*rand(3,nSpins));
                rnd = rnd(:,:)./sqrt(sum(rnd(:,:).*rnd(:,:)));

                dz = rnd*D;
                temp_r = dz+r;
                r = InBounds(constraint_radii(radius_bounds),temp_r(1,:),temp_r(2,:),temp_r(3,:), dz(1,:),dz(2,:),dz(3,:));
                Coords(radius_bounds,1,:, j) = r(1,:)';
                Coords(radius_bounds,2,:, j) = r(2,:)';
                Coords(radius_bounds,3,:, j) = r(3,:)';

                % condition is true if it has reached the read point
                if (j > diffusionGradient2_loc(end))
                    break;
                end

            end


    %     fig_num = 10;

    %     if (radius_bounds > 4)
    %         DrawParticleGraph(squeeze(xCoords(radius_bounds,:,:)),squeeze(yCoords(radius_bounds,:,:)),squeeze(zCoords(radius_bounds,:,:)),nSpins,fig_num)
    %     end

    end
save(sprintf('3D_Coords/Coords%d.mat',set),'Coords');
end
% save('xCoords20.mat','xCoords')
% save('yCoords20.mat','yCoords')
% save('zCoords20.mat','zCoords')


%% Perform sequence 20 times and continually add values

for SpinSet = 1:20
    
    fprintf("Set %d/20\n",SpinSet);
%     load(sprintf("Coords/xCoords%d.mat",SpinSet));
%     load(sprintf("Coords/yCoords%d.mat",SpinSet));
%     load(sprintf("Coords/zCoords%d.mat",SpinSet));
    load(sprintf("3D_Coords/Coords%d.mat",SpinSet));

    
    for radius_bounds = 1:length(constraint_radii)

        disp("Starting sequence");

        for i = 1:nG

            j = 1;
            gradAmp(3,diffusionGradient1_loc) =  G(i); %Z gradients in Tesla per meter
            gradAmp(3,diffusionGradient2_loc) =  G(i); %Z gradients in Tesla per meter

            mT = zeros(3,nSpins);
            mZ = ones(3,nSpins);

            %starting spin locations
            [mT,mZ] =  bloch(dt,([squeeze(Coords(radius_bounds,1,:,1)),...
                        squeeze(Coords(radius_bounds,2,:,1)),squeeze(Coords(radius_bounds,3,:,1))]'),0,T1,T2,mT,mZ); 

            for j = 2:nTimeSteps

                NOISE = 0;
                dB0 = gradAmp(:,j)'*([squeeze(Coords(radius_bounds,1,:,j)),...
                        squeeze(Coords(radius_bounds,2,:,j)),squeeze(Coords(radius_bounds,3,:,j))]') ...
                        + NOISE ; 

                [mT,mZ] =  bloch(dt,dB0,rfPulse(j),T1,T2,mT,mZ); 

                % condition is true if it has reached the read point
                if (j > diffusionGradient2_loc(end))
                    mFinalVect(i,:) = [mean(mT,'all'), mean(mZ,'all')];  
                    break;
                end

            end

            disp(['b-value = ' num2str(round(b(i)*1e-6))]);
        end


        store_final_vect(SpinSet,:,radius_bounds) = (mFinalVect(:,1));
    end
        
end
% legend("1*1e-6","10*1e-6","15*1e-6","20*1e-6","25*1e-6", "30*1e-6", "40*1e-6", "Free Space")
save('Large_final_vect.mat','store_final_vect')



%% Sum all results 
% save('Large_final_vect.mat','store_final_vect')
final_res = squeeze(sum((store_final_vect),1)./SetsOfSpins);
%%

for i = 1:length(constraint_radii)
    figure(2);
    subplot(2,1,1);
    log_vals = -log(abs(final_res(:,i))./(b*1e-6)');
    log_vals(1) = 1;
    plot(b*1e-6,log_vals,'.-','LineWidth',2);
    xlabel('b (s/mm^2)'), ylabel('ln(|M_{xy}|/b)');grid on;
    title('Log Diffusion Attenuation');
    hold on
    drawnow
    
    figure(2);
    subplot(2,1,2);
    plot(b*1e-6,abs(final_res(:,i)),'.-','LineWidth',2);
    xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
    title('Absolute of Diffusion Attenuation');
    hold on
    drawnow
end
figure(2);
subplot(2,1,1);
legend("1*1e-6","5*1e-6","10*1e-6","15*1e-6","20*1e-6","25*1e-6", "35*1e-6", "Free Space")
subplot(2,1,2);
legend("1*1e-6","5*1e-6","10*1e-6","15*1e-6","20*1e-6","25*1e-6", "35*1e-6", "Free Space")


%% FUNCTIONS %%

function result = InBounds(r,x1,y1,z1,dx,dy,dz)
    %result is dim = 3 x nSpins
    
    origin = [0 0 0];
    result = [x1;y1;z1];

    mag = sqrt((abs(x1)-origin(1)).^2 + (abs(y1)-origin(2)).^2 + (abs(z1)-origin(3)).^2);
    test_axis_x = sqrt((abs(x1)-origin(1)).^2 + (abs(y1-dy)-origin(2)).^2 + (abs(z1-dz)-origin(3)).^2);
    test_axis_y = sqrt((abs(x1-dx)-origin(1)).^2 + (abs(y1)-origin(2)).^2 + (abs(z1-dy)-origin(3)).^2);
    test_axis_z = sqrt((abs(x1-dx)-origin(1)).^2 + ((abs(y1-dy))-origin(2)).^2 + (abs(z1)-origin(3)).^2);
    test_axis = [test_axis_x;test_axis_y;test_axis_z];
    test_res = (test_axis>=r);
    temp = max(test_res);
    temp = max(temp);
    
    
    
    if (temp == 1)
        
        need_to_mirror = sum(test_res);
        index_list = find(need_to_mirror>0);
    
        mirrored = mirror_trajectory([(x1-dx);(y1-dy);(z1-dz)], [dx;dy;dz],index_list);
        result = mirrored;
    end
end


%TODO:
%Should probably vectorise this

function result = mirror_trajectory(previous_xyz,current_dxdydz, index_list)
    
    x0 = previous_xyz(1,:);
    y0 = previous_xyz(2,:);
    z0 = previous_xyz(3,:);
    
    x1 = x0 + current_dxdydz(1,:);
    y1 = y0 + current_dxdydz(2,:);
    z1 = z0 + current_dxdydz(3,:);
    
    result = [x1;y1;z1];
    
    vertex_x0y0z0 = [x0',y0',z0']; 
    vertex_x1y1z1 = [x1',y1',z1'];
    
%     figure;
    for i =1:length(index_list)
        
        V = [vertex_x0y0z0(index_list(i),:);...
             vertex_x1y1z1(index_list(i),:)];
         
%         plot3(V(:,1),V(:,2),V(:,3),'k.-','MarkerSize',25,'Color', rand(1,3)); 
%         xlabel('x'), ylabel('y'),zlabel('z');
%         hold on
        
        A = vertex_x0y0z0(index_list(i),:);
        B = vertex_x1y1z1(index_list(i),:);
        x = [V(1,1);V(2,1)];
        y = [V(1,2);V(2,2)];
        z = [V(1,3);V(2,3)];
        
        normal = ([mean(x),mean(y),mean(z)] + null(A-B)');
        normal = (normal(:,:)./sqrt(sum(normal(:,:).*normal(:,:))))*1e-6;%
        
        
%         plot3([B(1),normal(1)],[B(2),normal(2)],[B(3),normal(3)],'r.-','LineWidth',1)
%         legend("Original","Rotated")

        x1(index_list(i)) = normal(1);
        y1(index_list(i)) = normal(2);
        z1(index_list(i)) = normal(3);
        result = [x1;y1;z1];
    end

end

%Draws the plot for the particles trajetories 
function DrawParticleGraph(xCoords,yCoords,zCoords,nSpins,fig_num)
    figure(fig_num)
    for spin = 1:nSpins
        plot3(xCoords(spin,:),yCoords(spin,:),zCoords(spin,:),'Color', rand(1,3), 'MarkerSize', 9);
        hold on; 
    end
    xlabel('x'), ylabel('y'),zlabel('z');
    title('Final location of particle')
    drawnow;
end




function result = vector_mirror_trajectory(previous_xyz,current_dxdydz, index_list)
    
    x0 = previous_xyz(1,:);
    y0 = previous_xyz(2,:);
    z0 = previous_xyz(3,:);
    
    x1 = x0 + current_dxdydz(1,:);
    y1 = y0 + current_dxdydz(2,:);
    z1 = z0 + current_dxdydz(3,:);
    
    result = [x1;y1;z1];
    
    vertex_x0y0z0 = [x0',y0',z0']; 
    vertex_x1y1z1 = [x1',y1',z1'];
    
%     figure;
    for i =1:length(index_list)
        
        V = [vertex_x0y0z0(index_list(i),:);...
             vertex_x1y1z1(index_list(i),:)];
         
%         plot3(V(:,1),V(:,2),V(:,3),'k.-','MarkerSize',25,'Color', rand(1,3)); 
%         xlabel('x'), ylabel('y'),zlabel('z');
%         hold on
        
        A = vertex_x0y0z0(index_list(i),:);
        B = vertex_x1y1z1(index_list(i),:);
        x = [V(1,1);V(2,1)];
        y = [V(1,2);V(2,2)];
        z = [V(1,3);V(2,3)];
        
        normal = ([mean(x),mean(y),mean(z)] + null(A-B)');
        normal = (normal(:,:)./sqrt(sum(normal(:,:).*normal(:,:))))*1e-6;%
        
        
%         plot3([B(1),normal(1)],[B(2),normal(2)],[B(3),normal(3)],'r.-','LineWidth',1)
%         legend("Original","Rotated")

        x1(index_list(i)) = normal(1);
        y1(index_list(i)) = normal(2);
        z1(index_list(i)) = normal(3);
        result = [x1;y1;z1];
    end

end