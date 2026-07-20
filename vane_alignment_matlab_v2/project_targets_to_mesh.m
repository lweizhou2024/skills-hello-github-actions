function [measuredPoints,info] = project_targets_to_mesh( ...
        targetPoints,normals,TR,maxDistance,useNearestFallback)
%PROJECT_TARGETS_TO_MESH Intersect target-normal rays with a mesh.
%
% Rays are cast in both +normal and -normal directions. The nearest valid hit
% is selected. If neither ray hits, an optional closest-surface fallback is used.

    if nargin < 4 || isempty(maxDistance)
        maxDistance = Inf;
    end
    if nargin < 5 || isempty(useNearestFallback)
        useNearestFallback = true;
    end

    targetPoints = double(targetPoints);
    normals = normalizeRows(double(normals));

    V = TR.Points;
    F = TR.ConnectivityList;

    numberOfPoints = size(targetPoints,1);
    measuredPoints = zeros(numberOfPoints,3);
    info.mode = strings(numberOfPoints,1);
    info.distance = zeros(numberOfPoints,1);
    info.hitFace = nan(numberOfPoints,1);

    for i = 1:numberOfPoints
        [plusPoint,plusDistance,plusFace] = rayMeshNearestHit( ...
            targetPoints(i,:),normals(i,:),V,F,maxDistance);

        [minusPoint,minusDistance,minusFace] = rayMeshNearestHit( ...
            targetPoints(i,:),-normals(i,:),V,F,maxDistance);

        plusValid = all(isfinite(plusPoint));
        minusValid = all(isfinite(minusPoint));

        if plusValid || minusValid
            if plusValid && (~minusValid || plusDistance <= minusDistance)
                measuredPoints(i,:) = plusPoint;
                info.distance(i) = plusDistance;
                info.hitFace(i) = plusFace;
                info.mode(i) = "normal_ray_plus";
            else
                measuredPoints(i,:) = minusPoint;
                info.distance(i) = minusDistance;
                info.hitFace(i) = minusFace;
                info.mode(i) = "normal_ray_minus";
            end

        elseif useNearestFallback
            [nearestPoint,nearestFace,nearestDistance] = ...
                closest_points_on_mesh(targetPoints(i,:),TR,2);

            measuredPoints(i,:) = nearestPoint;
            info.distance(i) = nearestDistance;
            info.hitFace(i) = nearestFace;
            info.mode(i) = "nearest_surface_fallback";

        else
            error(['No surface intersection for target %d within ', ...
                   'maxDistance = %.6g.'],i,maxDistance);
        end
    end
end


function [hit,tHit,faceID] = rayMeshNearestHit(origin,direction,V,F,maxDistance)
% Vectorized Moller-Trumbore ray/triangle intersection.

    direction = direction/norm(direction);

    v0 = V(F(:,1),:);
    v1 = V(F(:,2),:);
    v2 = V(F(:,3),:);

    edge1 = v1-v0;
    edge2 = v2-v0;

    directionArray = repmat(direction,size(F,1),1);
    h = cross(directionArray,edge2,2);
    determinant = sum(edge1.*h,2);

    tolerance = 1e-12;
    valid = abs(determinant) > tolerance;

    inverseDeterminant = zeros(size(determinant));
    inverseDeterminant(valid) = 1./determinant(valid);

    s = origin-v0;
    u = inverseDeterminant.*sum(s.*h,2);

    q = cross(s,edge1,2);
    v = inverseDeterminant.*sum(directionArray.*q,2);
    rayDistance = inverseDeterminant.*sum(edge2.*q,2);

    barycentricTolerance = 1e-9;
    valid = valid & ...
            u >= -barycentricTolerance & ...
            v >= -barycentricTolerance & ...
            (u+v) <= 1+barycentricTolerance & ...
            rayDistance >= 0 & ...
            rayDistance <= maxDistance;

    candidates = find(valid);

    if isempty(candidates)
        hit = [NaN NaN NaN];
        tHit = Inf;
        faceID = NaN;
        return;
    end

    [tHit,localIndex] = min(rayDistance(candidates));
    faceID = candidates(localIndex);
    hit = origin+tHit*direction;
end


function X = normalizeRows(X)
    rowLength = vecnorm(X,2,2);
    if any(rowLength < eps)
        error("A zero-length normal vector was supplied.");
    end
    X = X./rowLength;
end
