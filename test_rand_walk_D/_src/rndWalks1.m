function coords = rndWalks1(D,nSpins,nTimeSteps,dt)
  %% allocate memory
  coords = zeros(3,nSpins,nTimeSteps);  
  %% usefull constrants
  DConst = sqrt(2*D*dt);  %dt is 1.0e-6
  
  coords(:,:,1) = zeros(3,nSpins);

    for ii = 2:nTimeSteps   
        % run steps
        coords(:,:,ii) = coords(:,:,ii-1) + DConst*randn(3,nSpins);
    end
%    plot_walks(coords,nSpins);
end
