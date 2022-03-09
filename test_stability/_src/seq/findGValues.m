function result = findGValues(b,bDelta,sDelta)

    gamma  = 2*pi*42.577e6; % rad/s/T
    result = sqrt(b/(gamma^2 * sDelta^2 * (bDelta - sDelta/3)));

end