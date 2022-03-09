function plt = plot_walks(Coords3D,nSpins)
    %% Plot of final distribution %%%
    figure(99);hold off 

    for particle = 1:nSpins
        x_ = squeeze(Coords3D(1,particle,:));
        y_ = squeeze(Coords3D(2,particle,:));
        z_ = squeeze(Coords3D(3,particle,:));

        plot3(x_, y_, z_, 'Color', rand(1,3), 'MarkerSize', 9);
        xlabel('x'), ylabel('y'),zlabel('z');
        hold on; 
    end
    hold off 
    
end