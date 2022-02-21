function result = getGVecMRISeq(G,gradAmp,rfPulse,T1,T2,dG1_loc,dG2_loc,Spin_locs,nSpins,dt)

    vec_GradAmp = zeros(length(G),3,500,size(gradAmp,2));
    
    idx = 1:length(G);
    vec_GradAmp(idx,3,:,dG1_loc) = G(idx)' .* ones(24,1,500,length(dG1_loc));
    vec_GradAmp(idx,3,:,dG2_loc) = G(idx)' .* ones(24,1,500,length(dG2_loc));
    
    mT = zeros(length(G),3,nSpins);
    mZ = ones(length(G), 3,nSpins);
    
    newSpin_locs = repmat(Spin_locs(:,:,:,:),1,1,1,24);
    newSpin_locs = permute(newSpin_locs,[4 1 2 3]);
   

    %starting spin locations
%     [mT,mZ] =  bloch(dt,([(Spin_locs(1,:,1)),...
%                           (Spin_locs(2,:,1)), ...
%                           (Spin_locs(3,:,1))]'),0,T1,T2,mT,mZ); 
    [mT,mZ] =  bloch(dt,squeeze(newSpin_locs(:,:,:,1)),0,T1,T2,mT,mZ); 

%     for j = 2:nTimeSteps
    j = 2;    
    while(1)
        

%         dB0 = gradAmp(:,j)'*([(Spin_locs(1,:,j)),...
%                               (Spin_locs(2,:,j)), ...
%                               (Spin_locs(3,:,j))]'); 

%         dB0 = vec_GradAmp(:,:,:,j)'*squeeze(newSpin_locs(:,:,:,j)); 
        
%         dB0 = pagemtimes(vec_GradAmp(:,:,1,j)',squeeze(newSpin_locs(:,:,:,j)));

         for idx = 1:length(G) %this dB0 should be 24x500
             dB0(idx,:) = vec_GradAmp(idx,:,:,j)*squeeze(newSpin_locs(idx,:,:,j));
         end
                          
        [mT,mZ] =  Copy_of_bloch(dt,dB0,rfPulse(j),T1,T2,mT,mZ); 

        % condition is true if it has reached the read point
        if (j > (dG2_loc(end)-1))
            result = [squeeze(mean(mean(mT,2),3)), squeeze(mean(mean(mZ,2),3))];  
            break;
        end
        j = j + 1;
    end

end