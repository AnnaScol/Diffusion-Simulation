% Finds the intesection point between a vector and a Sphere.
% r is the radius of the Sphere
% P is a point on the line
% U is the trajectory  vector P is on
function intersection = vectorSphereIntersection(radius, P, U)

   norm_traj = U./norm(U);

    a = dot(norm_traj,norm_traj);
    b = 2*dot(norm_traj,P);
    c = dot(P,P) - radius^2;
    d = b^2 - 4*a.*c;
    
    if d < 0
        fprintf("no intersection")
    else
        r = roots([a b c]);
        inters_magnitude = max(r);
        intersection = P + inters_magnitude*norm_traj;    

    end
    

end
