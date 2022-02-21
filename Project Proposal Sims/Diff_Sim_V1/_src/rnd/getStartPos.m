function pos = getStartPos(radius,nSpins)
 
    pos = 2*rand(3,nSpins*10)-1;
    
    r = dot(pos,pos,1)<1;
    
    pos = pos(:,r);
    
    pos = 0.999*radius*pos(:,1:nSpins);
    disp(size(pos));
            
%      figure; histogram(sqrt(dot(pos,pos,1)), 50); % should look exponential towwaars the egde

end