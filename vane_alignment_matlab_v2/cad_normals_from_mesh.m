function [normals,info] = cad_normals_from_mesh( ...
        targetPoints,TR,ringDepth,maxAngleDeg)
%CAD_NORMALS_FROM_MESH Estimate local CAD normals at target coordinates.
%
% The nearest CAD triangle supplies a seed normal. Nearby face normals are
% sign-aligned and area-averaged only when they are within maxAngleDeg of the
% seed. This avoids averaging across a sharp edge.

    if nargin < 3 || isempty(ringDepth)
        ringDepth = 1;
    end
    if nargin < 4 || isempty(maxAngleDeg)
        maxAngleDeg = 20;
    end

    [surfacePoint,seedFace,distanceToCAD] = ...
        closest_points_on_mesh(targetPoints,TR,max(1,ringDepth));

    V = TR.Points;
    F = TR.ConnectivityList;

    edge1 = V(F(:,2),:)-V(F(:,1),:);
    edge2 = V(F(:,3),:)-V(F(:,1),:);
    rawNormal = cross(edge1,edge2,2);
    doubleArea = vecnorm(rawNormal,2,2);
    faceNormal = rawNormal./max(doubleArea,eps);

    numberOfPoints = size(targetPoints,1);
    normals = zeros(numberOfPoints,3);
    angularSpreadDeg = zeros(numberOfPoints,1);
    numberOfFacesUsed = zeros(numberOfPoints,1);

    cosineLimit = cosd(maxAngleDeg);

    for i = 1:numberOfPoints
        seed = faceNormal(seedFace(i),:);
        seed = seed/norm(seed);

        activeVertices = unique(F(seedFace(i),:));
        candidateFaces = seedFace(i);

        for ring = 1:max(1,ringDepth)
            attached = vertexAttachments(TR,activeVertices);
            added = [];
            for k = 1:numel(attached)
                added = [added;attached{k}(:)]; %#ok<AGROW>
            end
            candidateFaces = unique([candidateFaces(:);added(:)]);
            activeVertices = unique(F(candidateFaces,:));
        end

        candidateNormals = faceNormal(candidateFaces,:);

        % STL triangle orientations can be inconsistent. Align every candidate
        % normal to the seed before averaging.
        directionSign = sign(candidateNormals*seed.');
        directionSign(directionSign == 0) = 1;
        candidateNormals = candidateNormals.*directionSign;

        cosineToSeed = candidateNormals*seed.';
        keep = cosineToSeed >= cosineLimit;

        candidateFaces = candidateFaces(keep);
        candidateNormals = candidateNormals(keep,:);

        if isempty(candidateFaces)
            localNormal = seed;
            spread = 0;
            count = 1;
        else
            weights = doubleArea(candidateFaces);
            localNormal = sum(candidateNormals.*weights,1);
            localNormal = localNormal/norm(localNormal);

            cosineToLocal = max(-1,min(1,candidateNormals*localNormal.'));
            spread = max(acosd(cosineToLocal));
            count = numel(candidateFaces);
        end

        normals(i,:) = localNormal;
        angularSpreadDeg(i) = spread;
        numberOfFacesUsed(i) = count;
    end

    info.surfacePoint = surfacePoint;
    info.seedFace = seedFace;
    info.distanceToCAD = distanceToCAD;
    info.angularSpreadDeg = angularSpreadDeg;
    info.numberOfFacesUsed = numberOfFacesUsed;
    info.ringDepth = ringDepth;
    info.maxAngleDeg = maxAngleDeg;
end
