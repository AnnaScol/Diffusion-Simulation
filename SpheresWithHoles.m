clear all; clc;

MainSphereRadius = 10;
MainSphereOrigin = [0,0,0];

CuttOutRadius = 2.5;
% CuttOutOrigin the origin will be the place along the mainsphere edge
N=100;
phi=linspace(0,pi,N);
theta=linspace(0,2*pi,2*N);
[phi,theta]=meshgrid(phi,theta);

cutout_disk = makeDisk(CuttOutRadius,N,N/2,2*N,N);
phi = phi.*(abs(1-cutout_disk));
theta = theta.*(abs(1-cutout_disk));

% Main sphere
x=MainSphereRadius.*sin(phi).*cos(theta);
y=MainSphereRadius.*sin(phi).*sin(theta);
z=MainSphereRadius.*cos(phi);

figure;
imshow(cutout_disk,[])
figure;
mesh(x,y,z)


% %% Make Disk
function disk_matrix = makeDisk(radiusPixels,x_center,y_center,sizeMatrix_y,sizeMatrix_x)
    disk_matrix = zeros(sizeMatrix_y,sizeMatrix_x); % create square matrix of zeroes
    origin = [y_center x_center]; % "center" of the matrix because we know it is sizeMatrix by sizeMatrix
    circle_centered = (1:sizeMatrix_y)-origin(1);
%     circle_centered_x = (1:sizeMatrix_x)-origin(1);

    [xx,yy] = meshgrid(circle_centered,circle_centered); % create x and y grid
    
    disk_matrix(sqrt(xx.^2 + yy.^2) <= radiusPixels) = 1; % set points inside the radius equal to one
end