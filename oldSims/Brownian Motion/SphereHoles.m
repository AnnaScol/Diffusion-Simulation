function [sphere_X,sphere_Y,sphere_Z,cutout_disk] = SphereHoles(MainSphereRadius,CutOutRadius,MainSphereOrigin,N)
    % CuttOutOrigin the origin will be the place along the mainsphere edge
    phi=linspace(0,pi,N);
    theta=linspace(0,2*pi,2*N);
    [phi,theta]=meshgrid(phi,theta);

    % %% Make Disk
    cutout_disk = zeros(2*N,N); % create square matrix of zeroes
    origin = [N/2 N]; % "center" of the matrix because we know it is sizeMatrix by sizeMatrix
    circle_centered = (1:(2*N))-origin(1);
    [xx,yy] = meshgrid(circle_centered,circle_centered); % create x and y grid
    cutout_disk(sqrt(xx.^2 + yy.^2) <= CutOutRadius) = 1; % set points inside the radius equal to one


    phi = phi.*(abs(1-cutout_disk));
    theta = theta.*(abs(1-cutout_disk));

    % Main sphere
    sphere_X=MainSphereRadius.*sin(phi).*cos(theta);
    sphere_Y=MainSphereRadius.*sin(phi).*sin(theta);
    sphere_Z=MainSphereRadius.*cos(phi);

%     figure;
%     imshow(cutout_disk,[])
%     figure;
%     mesh(sphere_X,sphere_Y,sphere_Z)

end