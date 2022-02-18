function result = simulateMRISequence(G,gradAmp,rfPulse,T1,T2,dG1_loc,dG2_loc,Spin_locs,nSpins,dt)

    gradAmp(3,dG1_loc) =  G; %Z gradients in Tesla per meter
    gradAmp(3,dG2_loc) =  G; %Z gradients in Tesla per meter

    mT = zeros(3,nSpins);
    mZ = ones(3,nSpins);

    %starting spin locations
%     [mT,mZ] =  bloch(dt,([(Spin_locs(1,:,1)),...
%                           (Spin_locs(2,:,1)), ...
%                           (Spin_locs(3,:,1))]'),0,T1,T2,mT,mZ); 
    [mT,mZ] =  bloch(dt,Spin_locs(:,:,1),0,T1,T2,mT,mZ); 

%     for j = 2:nTimeSteps
    j = 2;    
    while(1)
        

%         dB0 = gradAmp(:,j)'*([(Spin_locs(1,:,j)),...
%                               (Spin_locs(2,:,j)), ...
%                               (Spin_locs(3,:,j))]'); 

        dB0 = gradAmp(:,j)'*Spin_locs(:,:,j); 
                          
        [mT,mZ] =  bloch(dt,dB0,rfPulse(j),T1,T2,mT,mZ); 

        % condition is true if it has reached the read point
        if (j > (dG2_loc(end)-1))
            result = [mean(mT,'all'), mean(mZ,'all')];  
            break;
        end
        j = j + 1;
    end

end