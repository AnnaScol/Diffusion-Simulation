function coords = getRandomWalks(D,cellRadius,nSpins,nTimeSteps,dt)
  %% allocate memory
  
  coords = zeros(3,nSpins,nTimeSteps);
  fubar  = 0.0;
  
  %% usefull constrants
  sd__   = 1;                  %sqrt(STOP_TIME/numberOfSteps); 5.7735e-04;%
  nD     = 3;                     % We allways have 3D in this program
  DConst = sqrt(2*nD*D*dt)*sd__;  %tuned parameter to obtain the correct
  esc    = 0; % max number or tries

  coords(:,:,1) = getStartPos(cellRadius,nSpins);
 

    for ii = 2:nTimeSteps   

        % run steps

        coords(:,:,ii) = coords(:,:,ii-1) + randn(3,nSpins)*DConst;
        
        fubar = squeeze(dot(coords(:,:,ii),coords(:,:,ii),1));
        fu  = sqrt(fubar)>cellRadius;
        bar = sum(fu(:));
        esc = 0;
        
        % fix steps outside
        while bar>0 && esc<100
            coords(:,fu,ii) = coords(:,fu,ii-1) +randn(3,bar)*DConst;

            fubar = squeeze(dot(coords(:,:,ii),coords(:,:,ii),1));
            fu  = sqrt(fubar)>cellRadius;
            bar = sum(fu(:));
            
            esc = esc+1;
        end
        
        if max(sqrt(fubar))>cellRadius
            disp(['OOPS = # ' , num2str(bar)]);
        end
        
    end

    
         %% Plot of final distribution %%%
%             figure(99);hold off 
%     
%             for particle = 1:nSpins
%                 x_ = squeeze(coords(1,particle,:));
%                 y_ = squeeze(coords(2,particle,:));
%                 z_ = squeeze(coords(3,particle,:));
%                 
%                 plot3(x_, y_, z_, 'Color', rand(1,3), 'MarkerSize', 9);
%                 xlabel('x'), ylabel('y'),zlabel('z');
%                 hold on; 
%             end
%             hold off
    
end

