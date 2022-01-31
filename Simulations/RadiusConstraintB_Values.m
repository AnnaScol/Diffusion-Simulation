%% Spin Echo Line Scan
%do a 1d , get log plot
%map to a paper
%try different strating points
%e^{b*distanct traveled}
%finding the phase change
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 10^-6; 
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

nSpins = 1000;
G = ([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,100]*1e-3); %mT


nG = length(G);
mFinalVect = zeros(nG,2); %variable to hold the final magnetization calculated for each position

b = gamma^2*G.^2*sdelta^2*(ldelta-sdelta/3);

T1 = 1000*(10^-3);
T2 = 1000*(10^-3); 
TE = 65*(10^-3);

% Make size soma cell and brain cell sizes - vary 0.005-0.1mm
% radius_values = linspace(0.000005,0.0001,0.000001);


% constraint_radius = 5.0e-5;
% constraint_radius = 17.0e-6;
constraint_radii = ([1 5 10 25 50]*1e-6); 
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
gradAmp(3,diffusionGradient1_loc) =  G(3); %Z gradients in Tesla per meter

%RF Refocusing pules
pulsedurR = 0.002; % duration of the RF in s
rfStepsR = round(1:(pulsedurR/(dt)));
rfPulseR = apodize_sinc_rf(length(rfStepsR),3,pi,dt); %B1+ in Tesla
rfPulse(round(TE/2/dt) +length(rfStepsE)/2 + rfStepsR) = rfPulseR;

% diffusion pulse 2
diffusionGradient2_loc = round((1:(sdelta/dt)) + pulsedurE/dt + pulsedurR/dt + ldelta/dt + sdelta/dt);
gradAmp(3,diffusionGradient2_loc) =  G(3); %Z gradients in Tesla per meter

location = zeros(3,nTimeSteps);

%%% Plotting the data points %%%
xCoords = zeros(nSpins,nTimeSteps); %particle start loc is assume 0,0,0
yCoords = zeros(nSpins,nTimeSteps);
zCoords = zeros(nSpins,nTimeSteps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% PLOTTING %% %
figure
subplot(3,1,1); plot(time,rfPulse,'k-','LineWidth',2);title('RF Pulse'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;

subplot(3,1,2); plot(time,phase(rfPulse),'k-','LineWidth',2);title('RF Pulse Phase'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;

subplot(3,1,3); plot(time,gradAmp(3,:),'b-','LineWidth',2);title('Slice Select Gradient');
xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;

%% Get location matrix --> 3 x nSpins dimensions

for radius_bounds = 1:length(constraint_radii)
    for i = 1:nG
        j = 1;

        %starting spin locations
        r = zeros(3,nSpins);
        xCoords(i, j) = r(1);
        yCoords(i, j) = r(2);
        zCoords(i, j) = r(3);

        for j = 2:nTimeSteps   
            %%% Bounds checking %%%
            rnd = (-1 + 2*rand(3,nSpins));
            rnd = rnd(:,:)./sqrt(sum(rnd(:,:).*rnd(:,:)));

            dz = rnd*D;
            temp_r = dz+r;
            r = InBounds(constraint_radii(radius_bounds),temp_r(1,:),temp_r(2,:),temp_r(3,:), dz(1,:),dz(2,:),dz(3,:));
            xCoords(:, j) = r(1,:)';
            yCoords(:, j) = r(2,:)';
            zCoords(:, j) = r(3,:)';

            % condition is true if it has reached the read point
            if (j > diffusionGradient2_loc(end))
                break;
            end

        end
    end

    figure(10)
    for spin = 1:nSpins
        plot3(xCoords(spin,:),yCoords(spin,:),zCoords(spin,:),'Color', rand(1,3), 'MarkerSize', 9);
        hold on; 
    end
    xlabel('x'), ylabel('y'),zlabel('z');
    title('Final location of particle')
    drawnow;

    

    %% Perform sequence

    for i = 1:nG
        j = 1;
        gradAmp(3,diffusionGradient1_loc) =  G(i); %Z gradients in Tesla per meter
        gradAmp(3,diffusionGradient2_loc) =  G(i); %Z gradients in Tesla per meter

        mT = zeros(3,nSpins);
        mZ = ones(3,nSpins);

        %starting spin locations
        [mT,mZ] =  bloch(dt,([xCoords(:,1),yCoords(:,1),zCoords(:,1)]'),0,T1,T2,mT,mZ); 

        for j = 2:nTimeSteps

            dB0 = gradAmp(:,j)'*([xCoords(:,j),yCoords(:,j),zCoords(:,j)]'); 
            [mT,mZ] =  bloch(dt,dB0,rfPulse(j),T1,T2,mT,mZ); 

            % condition is true if it has reached the read point
            if (j > diffusionGradient2_loc(end))
                mFinalVect(i,:) = [mean(mT,'all'), mean(mZ,'all')];  
                break;
            end

        end

        disp(['b-value = ' num2str(round(b(i)*1e-6))]);
    end
    
    figure(2);
    log_vals = -log(abs(mFinalVect(:,1))./(b*1e-6)');
    log_vals(1) = 1;
    plot(b*1e-6,log_vals,'o-','LineWidth',2);
    xlabel('b (s/mm^2)'), ylabel('ln(|M_{xy}|/b)');grid on;
    title('Log Diffusion Attenuation');
    drawnow
    hold on;
    
end
legend("1*1e-6","5*1e-6","10*1e-6","25*1e-6","50*1e-6")
%0.02
%% PLOTTING %%
% Make it log values --> -(1/b)*ln(S_{DWI}/S_{b}) (dont need th 1/b maybe)
%make sure to scale both axis
%

figure;
log_vals = -log(abs(mFinalVect(:,1))./(b*1e-6)');
log_vals(1) = 1;
subplot(4,1,1);plot(b*1e-6,log_vals,'o-','LineWidth',2);
xlabel('b (s/mm^2)'), ylabel('ln(|M_{xy}|/b)');grid on;
title('Log Diffusion Attenuation');

subplot(4,1,2);plot(b*1e-6,abs(mFinalVect(:,1)),'o-','LineWidth',2);
xlabel('b (s/mm^2)'), ylabel('|M_{xy}|');grid on;
title('Absolute of Diffusion Attenuation');

% subplot(4,1,3);plot(b*1e-6,angle(mFinalVect(:,1)),'o-','LineWidth',2);
% xlabel('b (s/mm^2)'), ylabel('\angleM_{xy} '); ylim([-pi pi]);grid on;
% title('Angle of Diffusion Attenuation')

% 
% subplot(4,1,4);
figure
for spin = 1:nSpins
    plot3(xCoords(spin,:),yCoords(spin,:),zCoords(spin,:),'Color', rand(1,3), 'MarkerSize', 9);
    hold on; 
end
xlabel('x'), ylabel('y'),zlabel('z');
title('Final location of particle')



%% FUNCTIONS %%

function result = InBounds(r,x1,y1,z1,dx,dy,dz)
    %result is dim = 3 x nSpins
    
    origin = [0 0 0];
%     result = [x1;y1;z1];

    mag = sqrt((abs(x1)-origin(1)).^2 + (abs(y1)-origin(2)).^2 + (abs(z1)-origin(3)).^2);
    test_axis_x = sqrt((abs(x1)-origin(1)).^2 + (abs(y1-dy)-origin(2)).^2 + (abs(z1-dz)-origin(3)).^2);
    test_axis_y = sqrt((abs(x1-dx)-origin(1)).^2 + (abs(y1)-origin(2)).^2 + (abs(z1-dy)-origin(3)).^2);
    test_axis_z = sqrt((abs(x1-dx)-origin(1)).^2 + ((abs(y1-dy))-origin(2)).^2 + (abs(z1)-origin(3)).^2);
    test_axis = [test_axis_x;test_axis_y;test_axis_z];
    test_res = (test_axis>=r);
    temp = max(test_res);
    temp = max(temp);
    
    if (temp == 1)
        
        if (find(test_res(1,:)>0))
            x1(find(test_res(1,:)>0)) = x1(find(test_res(1,:)>0)) - dx(find(test_res(1,:)>0));
        end
        
        if (find(test_res(2,:)>0))
            y1(find(test_res(2,:)>0)) = y1(find(test_res(2,:)>0)) - dy(find(test_res(2,:)>0));
        end
        
        if (find(test_res(3,:)>0))
            z1(find(test_res(3,:)>0)) = z1(find(test_res(3,:)>0)) - dz(find(test_res(3,:)>0));
        end
    end
    
    result = [x1;y1;z1];
end
