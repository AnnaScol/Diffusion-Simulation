function [Sx,Sy,Sz] = CutOutSphere(MainSphereRadius,N)
%CuttOutRadius and CuttOutCenter are vectors of all the values
%N is number of separation points for phi and theta
% all areas of the cutouts are defined as NaN in the Sx,Sy,Sz

    % CuttOutOrigin the origin will be the place along the mainsphere edge
    phi=linspace(0,pi,N);
    theta=linspace(0,2*pi,N);
    [phi,theta]=meshgrid(phi,theta);


    % Make main sphere
    x=MainSphereRadius.*sin(phi).*cos(theta);
    y=MainSphereRadius.*sin(phi).*sin(theta);
    z=MainSphereRadius.*cos(phi);


    Sx = x; Sy = y; Sz = z;
    
    mesh(x,y,z);

end
