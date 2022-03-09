function coords = rndWalks2(D,nSpins,nTimeSteps,dt)
  %% allocate memory
  coords = zeros(3,nSpins,nTimeSteps);  
  %% usefull constrants
%   DConst = sqrt(2*D*dt);  %dt is 1.0e-6
  l = sqrt(2*D*3*0.01)% sqrt(2*D*0.001);
% predMean = sqrt(2*D*3*duration);
  coords(:,:,1) = zeros(3,nSpins);

    for ii = 2:nTimeSteps   
        % run steps
        step = 1-2*rand(3,nSpins);
        dxdydz = (step./vecnorm(step))*l;
        
        coords(:,:,ii) = coords(:,:,ii-1) + dxdydz;
        
    end
%    plot_walks(coords,nSpins);
end
