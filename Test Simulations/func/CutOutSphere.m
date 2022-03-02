function [Sx,Sy,Sz,cutout_disks,cutout_idx] = CutOutSphere(MainSphereRadius,CuttOutRadius,CuttOutCenter,N)
%CuttOutRadius and CuttOutCenter are vectors of all the values
%N is number of separation points for phi and theta
% all areas of the cutouts are defined as NaN in the Sx,Sy,Sz

    % CuttOutOrigin the origin will be the place along the mainsphere edge
    phi=linspace(0,pi,N);
    theta=linspace(0,2*pi,N);
    [phi,theta]=meshgrid(phi,theta);

    cutout_disks = zeros(length(CuttOutRadius),N,N);
    % cutout_idx = zeros(length(CuttOutRadius),200);
    for i = 1:length(CuttOutRadius)
        cutout_disks(i,:,:) = makeDisk(CuttOutRadius(i),N,abs(CuttOutCenter(i)),N,N);
        if(CuttOutCenter(i)< 0)
            cutout_disks(i,:,:) = flip(cutout_disks(i,:,:));
        end
        cutout_idx(i,1:length(find(squeeze(cutout_disks(i,:,:))==1))) = find(squeeze(cutout_disks(i,:,:))==1);

    end

    % Make main sphere
    x=MainSphereRadius.*sin(phi).*cos(theta);
    y=MainSphereRadius.*sin(phi).*sin(theta);
    z=MainSphereRadius.*cos(phi);


    Sx = x; Sy = y; Sz = z;
    
%     for i = 1:length(CuttOutRadius)
%         x(cutout_idx(i,(find(cutout_idx(i,:)>1)))) = nan;
%         y(cutout_idx(i,(find(cutout_idx(i,:)>1)))) = nan;
%         z(cutout_idx(i,(find(cutout_idx(i,:)>1)))) = nan; 
%     end
%     figure;
%     mesh(x,y,z);

end


% %% Make Disk
function disk_matrix = makeDisk(radiusPixels,x_center,y_center,sizeMatrix_y,sizeMatrix_x)
    disk_matrix = zeros(sizeMatrix_y,sizeMatrix_x); % create square matrix of zeroes
    origin = [y_center x_center]; % "center" of the matrix because we know it is sizeMatrix by sizeMatrix
    circle_centered = (1:sizeMatrix_y)-origin(1);
    [xx,yy] = meshgrid(circle_centered,circle_centered); % create x and y grid
    
    disk_matrix(sqrt(xx.^2 + yy.^2) <= radiusPixels) = 1; % set points inside the radius equal to one
end