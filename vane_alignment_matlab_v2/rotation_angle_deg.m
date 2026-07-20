function angleDeg = rotation_angle_deg(R)
%ROTATION_ANGLE_DEG Rotation magnitude represented by a 3-by-3 matrix.
    cosineValue = max(-1,min(1,(trace(R)-1)/2));
    angleDeg = acosd(cosineValue);
end
