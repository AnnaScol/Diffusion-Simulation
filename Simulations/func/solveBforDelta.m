function result = solveBforDelta(what_delta, other_delta_duration, b_values ,GradAmp)
    results = zeros(1,length(b_values));
    gamma = 2*pi*42.577*10^6;
    
    %find large delta for a specfic b-values and gradient
%     b_values = [1 25 50 100 200 300 400 500 750 1000 1500 2000 2500 3000 4000 ...
%                     5000 6000 7000 8000 9000 10000 11000 13000 14000 16500 18000]/1e-6;
                
    if (what_delta == 'b')
        for i = 1:length(b_values)
            partial_sol = b_values(i)/((gamma^2)*GradAmp^2);
            bDelta = (partial_sol/(other_delta_duration^2))+ (other_delta_duration/3);
            results(i) = bDelta;
        end
    elseif (what_delta == 'b')
        for i = 1:length(b_values)
            partial_sol = b_values(i)/((gamma^2)*GradAmp^2);
            bDelta = (partial_sol/(other_delta_duration^2))+ (other_delta_duration/3);
            results(i) = bDelta;
        end
    else
        fprintf("WRONG WHAT DELTA VALUE (must be 'b' for now)\n")
    end

%     disp(results);
    
    result = results;

end
