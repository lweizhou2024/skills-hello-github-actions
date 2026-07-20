function [TR, info] = normalize_stl_mesh(meshInput, mergeTolerance)
%NORMALIZE_STL_MESH Convert common STL/mesh representations to triangulation.
%
% [TR,INFO] = normalize_stl_mesh(meshInput, mergeTolerance)
%
% meshInput may be:
%   - MATLAB triangulation object
%   - filename accepted by stlread
%   - struct with any of these field pairs:
%       faces / vertices
%       face  / vertex
%       ConnectivityList / Points
%
% The function merges duplicated vertices. This is especially useful for STL
% structures in which each triangle owns three new vertex IDs even when the
% coordinates are identical to vertices in neighboring triangles.
%
% mergeTolerance:
%   0 or omitted : merge exactly equal coordinates
%   > 0          : merge coordinates after quantization by this tolerance
%
% Output:
%   TR   : MATLAB triangulation object
%   INFO : conversion and cleanup counts
%
% Example:
%   [TR,info] = normalize_stl_mesh(exchangeSTL,0);
%   clean.faces    = TR.ConnectivityList;
%   clean.vertices = TR.Points;

    if nargin < 2 || isempty(mergeTolerance)
        mergeTolerance = 0;
    end
    validateattributes(mergeTolerance,{'numeric'},{'scalar','nonnegative','finite'});

    sourceType = "";
    if isa(meshInput,'triangulation')
        F = meshInput.ConnectivityList;
        V = meshInput.Points;
        sourceType = "triangulation";

    elseif ischar(meshInput) || (isstring(meshInput) && isscalar(meshInput))
        sourceType = "file";
        filename = char(meshInput);

        try
            raw = stlread(filename);
            if isa(raw,'triangulation')
                F = raw.ConnectivityList;
                V = raw.Points;
            elseif isstruct(raw)
                [F,V] = extractFacesVertices(raw);
            else
                error("Unsupported one-output stlread result.");
            end
        catch firstError
            try
                [F,V] = stlread(filename);
            catch
                rethrow(firstError);
            end
        end

    elseif isstruct(meshInput)
        sourceType = "struct";
        [F,V] = extractFacesVertices(meshInput);

    else
        error(['meshInput must be a triangulation, STL filename, or a struct ', ...
               'containing faces and vertices.']);
    end

    F = double(F);
    V = double(V);

    if size(F,2) ~= 3
        error("Faces must be an N-by-3 array of triangle vertex IDs.");
    end
    if size(V,2) ~= 3
        error("Vertices must be an M-by-3 array of XYZ coordinates.");
    end
    if isempty(F) || isempty(V)
        error("Mesh faces and vertices cannot be empty.");
    end
    if any(~isfinite(F(:))) || any(~isfinite(V(:)))
        error("Mesh contains NaN or Inf.");
    end
    if any(abs(F(:)-round(F(:))) > 0)
        error("Face vertex IDs must be integers.");
    end
    F = round(F);

    if min(F(:)) == 0 && max(F(:)) <= size(V,1)-1
        warning("Detected zero-based face IDs. Converting them to MATLAB one-based IDs.");
        F = F + 1;
    end

    if min(F(:)) < 1 || max(F(:)) > size(V,1)
        error("Face IDs are outside the valid vertex range.");
    end

    originalVertexCount = size(V,1);
    originalFaceCount = size(F,1);

    % Merge duplicate vertices.
    if mergeTolerance == 0
        [Vmerged,~,vertexMap] = unique(V,'rows');
    else
        origin = min(V,[],1);
        key = round((V-origin)/mergeTolerance);
        [~,~,vertexMap] = unique(key,'rows');

        nMerged = max(vertexMap);
        Vmerged = zeros(nMerged,3);
        for column = 1:3
            Vmerged(:,column) = accumarray(vertexMap,V(:,column), ...
                [nMerged,1],@mean);
        end
    end

    Fmerged = reshape(vertexMap(F(:)),size(F));

    % Remove triangles collapsed by vertex merging.
    sortedFaces = sort(Fmerged,2);
    validDistinct = sortedFaces(:,1) ~= sortedFaces(:,2) & ...
                    sortedFaces(:,2) ~= sortedFaces(:,3);

    Fmerged = Fmerged(validDistinct,:);

    % Remove zero-area triangles.
    e1 = Vmerged(Fmerged(:,2),:) - Vmerged(Fmerged(:,1),:);
    e2 = Vmerged(Fmerged(:,3),:) - Vmerged(Fmerged(:,1),:);
    twiceArea = vecnorm(cross(e1,e2,2),2,2);

    scale = norm(max(Vmerged,[],1)-min(Vmerged,[],1));
    areaTolerance = max(eps(max(scale,1)^2)*100,realmin);
    validArea = twiceArea > areaTolerance;
    Fmerged = Fmerged(validArea,:);

    if isempty(Fmerged)
        error("No valid triangles remain after mesh cleanup.");
    end

    TR = triangulation(Fmerged,Vmerged);

    info.sourceType = sourceType;
    info.originalVertexCount = originalVertexCount;
    info.mergedVertexCount = size(Vmerged,1);
    info.verticesRemoved = originalVertexCount-size(Vmerged,1);
    info.originalFaceCount = originalFaceCount;
    info.cleanedFaceCount = size(Fmerged,1);
    info.facesRemoved = originalFaceCount-size(Fmerged,1);
    info.mergeTolerance = mergeTolerance;

    fprintf("Mesh conversion: %d -> %d vertices; %d -> %d faces.\n", ...
        originalVertexCount,size(Vmerged,1), ...
        originalFaceCount,size(Fmerged,1));
end


function [F,V] = extractFacesVertices(s)
    if isfieldCaseInsensitive(s,'TR')
        candidate = getFieldCaseInsensitive(s,'TR');
        if isa(candidate,'triangulation')
            F = candidate.ConnectivityList;
            V = candidate.Points;
            return;
        end
    end

    faceNames = {'ConnectivityList','faces','face'};
    vertexNames = {'Points','vertices','vertex'};

    F = [];
    V = [];

    for i = 1:numel(faceNames)
        if isfieldCaseInsensitive(s,faceNames{i})
            F = getFieldCaseInsensitive(s,faceNames{i});
            break;
        end
    end

    for i = 1:numel(vertexNames)
        if isfieldCaseInsensitive(s,vertexNames{i})
            V = getFieldCaseInsensitive(s,vertexNames{i});
            break;
        end
    end

    if isempty(F) || isempty(V)
        error(['Mesh struct must contain faces/vertices, face/vertex, or ', ...
               'ConnectivityList/Points fields.']);
    end
end


function tf = isfieldCaseInsensitive(s,name)
    fields = fieldnames(s);
    tf = any(strcmpi(fields,name));
end


function value = getFieldCaseInsensitive(s,name)
    fields = fieldnames(s);
    index = find(strcmpi(fields,name),1);
    value = s.(fields{index});
end
