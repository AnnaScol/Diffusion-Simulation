%% Spin Echo Line Scan
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 10^-5; 
gamma = 2*pi*42.577*10^6;


%Allocate the memory needed
nTimeSteps  = 6000;%70ms
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
adc         = zeros(1,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points

% %% Parameters %% %

% Diffusion Gradients
ldelta = 0.030; %ms
sdelta = 0.015; %ms
D = 1.5e-12;    %m^2/s

nSpins = 500;
start_position = zeros(1,nSpins);
% G = ([0,5,10]*1e-3); %mT

G = ([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,110]*1e-3); %mT

nG = length(G);
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position

b = gamma^2*G.^2*sdelta^2*(ldelta-sdelta/3);

T1 = 850*(10^-3);
T2 = 80*(10^-3); 
TE = 60*(10^-3);

constraint_radius = 2.0e-5;
% constraint_radius = 5;

% Sphere center
x_center = 0; y_center = 0; z_center = 0;

% Generate Sphere with hole
N=100;  MainSphereOrigin = [x_center,y_center,z_center];
MainSphereRadius = constraint_radius;
CuttOutRadius = [5];
CuttOutCenter = [30];

[sphere_X,sphere_Y,sphere_Z,cutout_disk,cutout_idx] = CutOutSphere(MainSphereRadius,CuttOutRadius,CuttOutCenter,N);

% This gets passed into the boundary checker
boundary.x = sphere_X;
boundary.y = sphere_Y;
boundary.z = sphere_Z;
boundary.cutouts = cutout_disk;
boundary.cutout_idx = cutout_idx;
%%

% %%   %% %
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

% Generate and Display MRI Sequence
% %% RF %% %
% RF Excitation Waveform
pulsedurE = 0.001; % duration of the RF in s
rfStepsE = round(1:(pulsedurE/(dt)));
rfPulse(rfStepsE) = apodize_sinc_rf(length(rfStepsE),3,pi/2,dt); %B1+ in Tesla

% First diffusion gradient
diffusionGradient1_loc = round((1:(sdelta/dt))+ pulsedurE/dt);
gradAmp(3,diffusionGradient1_loc) =  G(2); %Z gradients in Tesla per meter

%RF Refocusing pules
pulsedurR = 0.001; % duration of the RF in s
rfStepsR = round(1:(pulsedurR/(dt)));
rfPulseR = apodize_sinc_rf(length(rfStepsR),3,pi,dt); %B1+ in Tesla
rfPulse(round(TE/2/dt) +length(rfStepsE)/2 + rfStepsR) = rfPulseR;

% diffusion pulse 2
diffusionGradient2_loc = round((1:(sdelta/dt)) + pulsedurE/dt + pulsedurR/dt + ldelta/dt);
gradAmp(3,diffusionGradient2_loc) =  G(2); %Z gradients in Tesla per meter


% %% PLOTTING %% %
figure
subplot(3,1,1); plot(time,rfPulse,'k-','LineWidth',2);title('RF Pulse'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;

subplot(3,1,2); plot(time,phase(rfPulse),'k-','LineWidth',2);title('RF Pulse Phase'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;

subplot(3,1,3); plot(time,gradAmp(3,:),'b-','LineWidth',2);title('Slice Select Gradient');
xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;

%% Perform seuquence

for i = 1:nG
    gradAmp(3,diffusionGradient1_loc) =  G(i); %Z gradients in Tesla per meter
    gradAmp(3,diffusionGradient2_loc) =  G(i); %Z gradients in Tesla per meter
    
    mT = zeros(1,nSpins);
    mZ = ones(1,nSpins);
    
    %starting spin locations
    deltaZ = zeros(3,nSpins);
    [mT,mZ] =  bloch(dt,deltaZ,0,T1,T2,mT,mZ); 
    
    for j = 2:nTimeSteps
        
        %%% Bounds checking %%%
        dz = (randn(3,nSpins)).*sqrt(2*D*(ldelta-sdelta/3));
        temp = dz+deltaZ; %adding a random step to each spin in each timestep between diffusion grads
        deltaZ = InBounds(constraint_radius,temp(1,:),temp(2,:),temp(3,:), dz(1,:),dz(2,:),dz(3,:),boundary);
        %%%%%%%%%%%%%%%%%%%%%%%
        
        dB0 = gradAmp(:,j).*deltaZ; 
        [mT,mZ] =  bloch(dt,dB0,rfPulse(j),T1,T2,mT,mZ); 
        
        % condition is true if it has reached the read point
        if (j > diffusionGradient2_loc(end))
            mFinalVect(i,:) = [mean(mT,'all'), mean(mZ,'all')];  
            break;
        end

    end
   
    disp(['b-value = ' num2str(round(b(i)*1e-6))]);
    
end


figure;
subplot(3,1,1);plot(b*1e-6,abs(mFinalVect(:,1)),'o-','LineWidth',2);
xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
title('Absolute of Diffusion Attenuation');

subplot(3,1,2);plot(b*1e-6,angle(mFinalVect(:,1)),'o-','LineWidth',2);
xlabel('b (s/mm^2)'), ylabel('\angleM_{xy} '); ylim([-pi pi]);grid on;
title('Angle of Diffusion Attenuation')


subplot(3,1,3);
for spin = 1:nSpins
    plot3(deltaZ(1,spin),deltaZ(2,spin),deltaZ(3,spin),'x');
    hold on; 
end
xlabel('x'), ylabel('y'),zlabel('z');
%plot line axis
line_ = -0.00005:0.000001:0.00005;
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

title('Final location of particle')


%% FUNCTIONS %%

function result = InBounds(r,x1,y1,z1,dx,dy,dz,sphere_info)
    checker = zeros(1,length(x1));
    origin = [0 0 0];
    mag = sqrt((abs(x1)-origin(1)).^2 + (abs(y1)-origin(2)).^2 + (abs(z1)-origin(3)).^2);
    outside_radius = (mag >= r);
    check_x1 = outside_radius.*x1;
    check_y1 = outside_radius.*y1;
    check_z1 = outside_radius.*z1;
    
        
   
    if (sum(check_x1)+sum(check_y1)+sum(check_z1) ~= 0)
        checker = CheckCutOut(origin,check_x1,check_y1,check_z1,sphere_info);
        
            if ((sum(checker) == 0) && (length(checker) == 1) && ~size(sphere_info.cutout_idx,2))
                
                index = find(outside_radius>0);
                x1(index) = x1(index)-dx(index); 
                y1(index) = y1(index)-dy(index); 
                z1(index) = z1(index)-dz(index); 
            end
            
            if (sum(checker)>0)
                index = find(checker>0);
                x1(index) = x1(index)-dx(index); 
                y1(index) = y1(index)-dy(index); 
                z1(index) = z1(index)-dz(index); 
            end
               
    end

    result = [x1;y1;z1];
end



function result = CheckCutOut(origin,x,y,z,sphere_info)

    sphere_X = sphere_info.x;
    sphere_Y = sphere_info.y;
    sphere_Z = sphere_info.z;
    cutout_disk = sphere_info.cutouts;
    cutout_idx = sphere_info.cutout_idx;
    result = zeros(1,length(x));

    %try testing each of the x,y,x and x,y,z sphere values in a loop
    for i = 1:size(cutout_idx,1)
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
            
            %maybe cascade the results downwards because it will only
            %correct one spins 
            
            %check if it is meant to be negative or positive in the axis limits
            if ( (Sx >= 0) && (Sy >= 0) && (Sz >= 0) )      
                if ( sum((Sx <= x) & (Sy <= y) & (Sz <= z)) > 0)
                    
                    result = ((Sx <= x) & (Sy <= y) & (Sz <= z)) ;
                end
                
            if ( (Sx >= 0) && (Sy >= 0) && (Sz < 0) )
                if ( (Sx <= x) & (Sy <= y) & (Sz >= z) )
                    
                    result = (Sx <= x) & (Sy <= y) & (Sz >= z);
                end
            end

            if ( (Sx >= 0) && (Sy < 0) && (Sz >= 0) )
                if ( (Sx <= x) & (Sy >= y) & (Sz <= z) )
                    
                    result = (Sx <= x) & (Sy >= y) & (Sz <= z) ;
                end
            end
                
            if ( (Sx < 0) && (Sy >= 0) && (Sz >= 0) )
                if ( (Sx >= x) & (Sy <= y) & (Sz <= z) )
                    
                    result = (Sx >= x) & (Sy <= y) & (Sz <= z) ;
                end 
            end
               
            if ( (Sx >= 0) && (Sy < 0) && (Sz < 0) )
                if ( (Sx <= x) & (Sy >= y) & (Sz >= z) )
                    
                    result = (Sx <= x) & (Sy >= y) & (Sz >= z) ;
                end  
            end
                
            if ( (Sx < 0) && (Sy < 0) && (Sz >= 0) )
                if ( (Sx  >= x) & (Sy >= y) & (Sz <= z) )
                    
                    result = (Sx  >= x) & (Sy >= y) & (Sz <= z) ;
                end 
            end
                
            if ( (Sx < 0) && (Sy >= 0) && (Sz < 0) )
                if ( (Sx  >= x) & (Sy <= y) & (Sz >= z) )
                    
                    result = (Sx  >= x) & (Sy <= y) & (Sz >= z) ;
                end 
            end
            
            if ( (Sx < 0) && (Sy < 0) && (Sz < 0) )
                if ( (Sx >= x) & (Sy >= y) & (Sz >= z) )
                    
                    result = (Sx >= x) & (Sy >= y) & (Sz >= z) ;
                end  
            end
                
            end %end of if,elseif,..... statement
        end %end of all index checks
    end %end of num_counts 

end



