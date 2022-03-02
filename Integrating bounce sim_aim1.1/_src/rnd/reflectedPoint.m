function reflected_coord = reflectedPoint(intersection_point, particle_location)

    reflected_coord = particle_location - 2*dot(particle_location,intersection_point)*intersection_point;

end
