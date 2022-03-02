function coords = rndWalks(D,nSpins,nTimeSteps,dt)
  %% allocate memory
  coords = zeros(3,nSpins,nTimeSteps);  
  %% usefull constrants
  nD     = 3;                     % We allways have 3D in this program
  DConst = sqrt(2*D*dt);  %tuned parameter to obtain the correct
  
  coords(:,:,1) = zeros(3,nSpins);

    for ii = 2:nTimeSteps   
        % run steps
        coords(:,:,ii) = coords(:,:,ii-1) + DConst*randn(3,nSpins);
    end
%    plot_walks(coords,nSpins);
end
