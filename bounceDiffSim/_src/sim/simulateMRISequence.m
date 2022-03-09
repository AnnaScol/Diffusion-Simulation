function result = simulateMRISequence(Sequence_Struct,T1,T2,spinCor,dt)

    nTimeSteps = length(Sequence_Struct.RF);
    nSpins     = size(spinCor,2);
        
    mT = zeros(1,nSpins);
    mZ = ones(1,nSpins);


    dB0 = Sequence_Struct.G(:,1)'*spinCor(:,:,1);
    [mT,mZ] =  bloch(dt,dB0,Sequence_Struct.RF(1),T1,T2,mT,mZ); 

     for ii = 2:nTimeSteps

        dB0 = Sequence_Struct.G(:,ii)'*spinCor(:,:,ii);                    
        [mT,mZ] =  bloch(dt,dB0,Sequence_Struct.RF(ii),T1,T2,mT,mZ);   
        
     end
     
     result = mean(mT);  

end