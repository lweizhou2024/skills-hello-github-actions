function [Q,faceID,distance] = closest_points_on_mesh(P,TR,ringDepth)
%CLOSEST_POINTS_ON_MESH Closest surface points using local triangle searches.
%
% The nearest mesh vertex is found first. Triangles in its local topological
% neighborhood are then tested exactly. ringDepth=2 is a good default.

    if nargin < 3 || isempty(ringDepth)
        ringDepth = 2;
    end

    P = double(P);
    if size(P,2) ~= 3
        error("P must be an N-by-3 coordinate array.");
    end

    V = TR.Points;
    F = TR.ConnectivityList;

    seedVertex = nearestNeighbor(TR,P);

    numberOfPoints = size(P,1);
    Q = zeros(numberOfPoints,3);
    faceID = zeros(numberOfPoints,1);
    distance = zeros(numberOfPoints,1);

    for i = 1:numberOfPoints
        activeVertices = seedVertex(i);
        candidateFaces = [];

        for ring = 1:max(1,ringDepth)
            attached = vertexAttachments(TR,activeVertices);
            newFaces = flattenCellNumbers(attached);
            candidateFaces = unique([candidateFaces(:);newFaces(:)]);

            if isempty(candidateFaces)
                break;
            end
            activeVertices = unique(F(candidateFaces,:));
        end

        if isempty(candidateFaces)
            error("Could not find triangles attached to nearest mesh vertex.");
        end

        bestDistanceSquared = Inf;
        bestPoint = [NaN NaN NaN];
        bestFace = NaN;

        for j = 1:numel(candidateFaces)
            currentFace = candidateFaces(j);
            triangle = V(F(currentFace,:),:);
            candidatePoint = closestPointTriangle( ...
                P(i,:),triangle(1,:),triangle(2,:),triangle(3,:));

            d2 = sum((candidatePoint-P(i,:)).^2);
            if d2 < bestDistanceSquared
                bestDistanceSquared = d2;
                bestPoint = candidatePoint;
                bestFace = currentFace;
            end
        end

        Q(i,:) = bestPoint;
        faceID(i) = bestFace;
        distance(i) = sqrt(bestDistanceSquared);
    end
end


function values = flattenCellNumbers(c)
    if isempty(c)
        values = [];
        return;
    end

    values = [];
    for i = 1:numel(c)
        values = [values;c{i}(:)]; %#ok<AGROW>
    end
    values = unique(values);
end


function q = closestPointTriangle(p,a,b,c)
% Ericson, Real-Time Collision Detection, closest point on triangle.

    ab = b-a;
    ac = c-a;
    ap = p-a;

    d1 = dot(ab,ap);
    d2 = dot(ac,ap);
    if d1 <= 0 && d2 <= 0
        q = a;
        return;
    end

    bp = p-b;
    d3 = dot(ab,bp);
    d4 = dot(ac,bp);
    if d3 >= 0 && d4 <= d3
        q = b;
        return;
    end

    vc = d1*d4-d3*d2;
    if vc <= 0 && d1 >= 0 && d3 <= 0
        v = d1/(d1-d3);
        q = a+v*ab;
        return;
    end

    cp = p-c;
    d5 = dot(ab,cp);
    d6 = dot(ac,cp);
    if d6 >= 0 && d5 <= d6
        q = c;
        return;
    end

    vb = d5*d2-d1*d6;
    if vb <= 0 && d2 >= 0 && d6 <= 0
        w = d2/(d2-d6);
        q = a+w*ac;
        return;
    end

    va = d3*d6-d5*d4;
    if va <= 0 && (d4-d3) >= 0 && (d5-d6) >= 0
        w = (d4-d3)/((d4-d3)+(d5-d6));
        q = b+w*(c-b);
        return;
    end

    denominator = 1/(va+vb+vc);
    v = vb*denominator;
    w = vc*denominator;
    q = a+v*ab+w*ac;
end
