function [ mtOut_nw, mzOut  ] = Copy_of_bloch(dt,dB0,B1,T1,T2,mt,mz) %1sec

   gamma = 2*pi*42.577*10^6; % rad Hz per Tesla
    
   %calulate mtOut and mzOut:
   mtOut = 1i*sin(gamma*dt*B1).*mz;
   mtOut = mtOut + cos(gamma*dt*imag(B1)).*real(mt) +1i*cos(gamma*dt*real(B1)).*imag(mt);
   
   mzOut = cos(gamma*B1*dt).*mz;
   mzOut = mzOut + sin(gamma*dt*imag(B1)).*real(mt) - sin(gamma*dt*real(B1)).*imag(mt);
   
   %apply B0:
   mtOut_nw = zeros(size(dB0,1),size(mtOut,1),size(mtOut,2));
   for idx = 1:size(dB0,1)
    mtOut_nw(idx,:,:) = mtOut.*exp(1i*gamma*dt*dB0(idx,:));
   end
   
%    %apply relaxation: 
%    mtOut = mtOut.*exp(-dt./T2);
%    mzOut = mzOut.*exp(-dt./T1) + (1-exp(-dt./T1));
end