function result = vane_rps_alignment(cfg)
%VANE_RPS_ALIGNMENT Datum-target alignment for nominal and distorted STL meshes.
%
% PRE-ALIGNED INPUT (recommended for your current workflow)
% ---------------------------------------------------------
% cfg.cadMesh       = nominalStruct;      % or triangulation / STL filename
% cfg.scanMesh      = prealignedStruct;   % or triangulation / STL filename
% cfg.preAligned    = true;
% cfg.datumPointsCAD = N-by-3 target coordinates from WENZEL;
% cfg.method        = "rps";              % or "point"
%
% OPTIONAL
% --------
% cfg.datumNormalsCAD       N-by-3. If omitted, derive from nominal CAD mesh.
% cfg.weights               N-by-1, default ones
% cfg.meshMergeTolerance    default 0, merges exact duplicate coordinates
% cfg.normalRingDepth       default 1
% cfg.normalMaxAngleDeg     default 20
% cfg.maxProjectionDistance default 5 model units
% cfg.maxIterations         default 8
% cfg.translationTolerance  default 1e-6
% cfg.rotationToleranceDeg  default 1e-6
% cfg.useNearestFallback    default true
% cfg.plot                  default true
%
% COARSE 8-CORNER ALIGNMENT (optional)
% ------------------------------------
% Set cfg.preAligned=false and supply:
% cfg.cadCorners
% cfg.scanCorners
%
% MESH INPUT
% ----------
% The mesh can use faces/vertices, face/vertex, ConnectivityList/Points, a
% triangulation object, or an STL filename. Repeated STL vertices are merged.
%
% Transformation convention:
% Xaligned = (result.R * Xoriginal.' + result.t).'

    cfg = defaults(cfg);

    cadInput = meshInputFromCfg(cfg,'cad');
    scanInput = meshInputFromCfg(cfg,'scan');

    [TRcad,cadMeshInfo] = normalize_stl_mesh( ...
        cadInput,cfg.meshMergeTolerance);
    [TRscan,scanMeshInfo] = normalize_stl_mesh( ...
        scanInput,cfg.meshMergeTolerance);

    targetPoints = double(cfg.datumPointsCAD);
    numberOfTargets = size(targetPoints,1);

    if size(targetPoints,2) ~= 3
        error("cfg.datumPointsCAD must be an N-by-3 array.");
    end

    weights = double(cfg.weights(:));
    if isempty(weights)
        weights = ones(numberOfTargets,1);
    end
    if numel(weights) ~= numberOfTargets || any(weights <= 0)
        error("cfg.weights must contain one positive value per target.");
    end
    weights = weights/mean(weights);

    if cfg.preAligned
        initialR = eye(3);
        initialT = zeros(3,1);
        initialCornerRMS = NaN;
    else
        requireField(cfg,'cadCorners');
        requireField(cfg,'scanCorners');
        [initialR,initialT,initialCornerRMS] = weighted_kabsch( ...
            cfg.scanCorners,cfg.cadCorners,[]);
    end

    if isfield(cfg,'datumNormalsCAD') && ~isempty(cfg.datumNormalsCAD)
        normals = normalizeRows(double(cfg.datumNormalsCAD));
        normalInfo.source = "supplied";
        normalInfo.distanceToCAD = nan(numberOfTargets,1);
        normalInfo.angularSpreadDeg = nan(numberOfTargets,1);
    else
        [normals,normalInfo] = cad_normals_from_mesh( ...
            targetPoints,TRcad,cfg.normalRingDepth,cfg.normalMaxAngleDeg);
        normalInfo.source = "derived_from_CAD_mesh";
    end

    if ~isequal(size(normals),size(targetPoints))
        error("Datum normals must match the size of datumPointsCAD.");
    end

    R = initialR;
    t = initialT;

    history = repmat(struct( ...
        'iteration',[], ...
        'deltaTranslation',[], ...
        'deltaRotationDeg',[], ...
        'normalRMS',[], ...
        'pointRMS',[]),cfg.maxIterations,1);

    for iteration = 1:cfg.maxIterations
        alignedVertices = apply_rigid(TRscan.Points,R,t);
        alignedTR = triangulation(TRscan.ConnectivityList,alignedVertices);

        [measuredPoints,projectionInfo] = project_targets_to_mesh( ...
            targetPoints,normals,alignedTR, ...
            cfg.maxProjectionDistance,cfg.useNearestFallback);

        if strcmpi(cfg.method,'rps')
            [incrementR,incrementT,fitInfo] = fit_rps_increment( ...
                measuredPoints,targetPoints,normals,weights);
        elseif any(strcmpi(cfg.method,{'point','point_to_point','kabsch'}))
            [incrementR,incrementT,pointRMS] = weighted_kabsch( ...
                measuredPoints,targetPoints,weights);

            fittedPoints = apply_rigid(measuredPoints,incrementR,incrementT);
            normalResidual = sum((fittedPoints-targetPoints).*normals,2);
            fitInfo.normalRMS = sqrt(sum(weights.*normalResidual.^2)/sum(weights));
            fitInfo.rank = NaN;
            fitInfo.conditionNumber = NaN;
        else
            error("cfg.method must be 'rps' or 'point'.");
        end

        R = incrementR*R;
        t = incrementR*t+incrementT;

        history(iteration).iteration = iteration;
        history(iteration).deltaTranslation = norm(incrementT);
        history(iteration).deltaRotationDeg = rotation_angle_deg(incrementR);
        history(iteration).normalRMS = fitInfo.normalRMS;

        fittedPoints = apply_rigid(measuredPoints,incrementR,incrementT);
        history(iteration).pointRMS = sqrt(sum(weights.* ...
            sum((fittedPoints-targetPoints).^2,2))/sum(weights));

        if history(iteration).deltaTranslation < cfg.translationTolerance && ...
                history(iteration).deltaRotationDeg < cfg.rotationToleranceDeg
            history = history(1:iteration);
            break;
        end
    end

    alignedVertices = apply_rigid(TRscan.Points,R,t);
    alignedTR = triangulation(TRscan.ConnectivityList,alignedVertices);

    [finalMeasuredPoints,projectionInfo] = project_targets_to_mesh( ...
        targetPoints,normals,alignedTR, ...
        cfg.maxProjectionDistance,cfg.useNearestFallback);

    normalResidual = sum((finalMeasuredPoints-targetPoints).*normals,2);
    pointResidual = vecnorm(finalMeasuredPoints-targetPoints,2,2);

    J = [cross(finalMeasuredPoints,normals,2),normals];
    weightedJ = sqrt(weights).*J;
    singularValues = svd(weightedJ,0);
    rankTolerance = max(size(weightedJ))*eps(max(singularValues));
    rpsRank = sum(singularValues > rankTolerance);

    if numel(singularValues) >= 6 && singularValues(end) > 0
        conditionNumber = singularValues(1)/singularValues(end);
    else
        conditionNumber = Inf;
    end

    result.R = R;
    result.t = t;
    result.initial.R = initialR;
    result.initial.t = initialT;
    result.initial.cornerRMS = initialCornerRMS;

    result.datumPointsCAD = targetPoints;
    result.datumNormalsCAD = normals;
    result.measuredDatumPointsAligned = finalMeasuredPoints;

    result.alignedScanTR = alignedTR;
    result.alignedScanMesh.faces = alignedTR.ConnectivityList;
    result.alignedScanMesh.vertices = alignedTR.Points;

    result.rps.normalResiduals = normalResidual;
    result.rps.pointResiduals = pointResidual;
    result.rps.normalRMS = sqrt(sum(weights.*normalResidual.^2)/sum(weights));
    result.rps.pointRMS = sqrt(sum(weights.*pointResidual.^2)/sum(weights));
    result.rps.rank = rpsRank;
    result.rps.singularValues = singularValues;
    result.rps.conditionNumber = conditionNumber;

    result.history = history;
    result.normalInfo = normalInfo;
    result.projectionInfo = projectionInfo;
    result.cadMeshInfo = cadMeshInfo;
    result.scanMeshInfo = scanMeshInfo;

    fprintf("\nVane datum alignment\n");
    fprintf("  Method                 : %s\n",string(cfg.method));
    fprintf("  Pre-aligned input      : %d\n",cfg.preAligned);
    fprintf("  Datum normal RMS       : %.6g\n",result.rps.normalRMS);
    fprintf("  Datum 3-D point RMS    : %.6g\n",result.rps.pointRMS);
    fprintf("  RPS constraint rank    : %d of 6\n",result.rps.rank);
    fprintf("  RPS condition number   : %.6g\n",result.rps.conditionNumber);
    fprintf("  Added rotation         : %.6g degrees\n",rotation_angle_deg(R));
    fprintf("  Added translation      : %.6g\n",norm(t));

    if result.rps.rank < 6
        warning(['Targets and normal directions do not constrain all six ', ...
                 'rigid-body degrees of freedom.']);
    elseif result.rps.conditionNumber > 1e6
        warning(['The RPS geometry is poorly conditioned. Review point spacing ', ...
                 'and normal directions.']);
    end

    if cfg.plot
        plotDatumAlignment(TRcad,alignedTR,targetPoints, ...
            finalMeasuredPoints,normals,cfg.method);
    end
end


function cfg = defaults(cfg)
    values.preAligned = [];
    values.method = "rps";
    values.weights = [];
    values.meshMergeTolerance = 0;
    values.normalRingDepth = 1;
    values.normalMaxAngleDeg = 20;
    values.maxProjectionDistance = 5;
    values.maxIterations = 8;
    values.translationTolerance = 1e-6;
    values.rotationToleranceDeg = 1e-6;
    values.useNearestFallback = true;
    values.plot = true;

    names = fieldnames(values);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(cfg,name) || isempty(cfg.(name))
            cfg.(name) = values.(name);
        end
    end

    if isempty(cfg.preAligned)
        cfg.preAligned = ~(isfield(cfg,'cadCorners') && ...
                           isfield(cfg,'scanCorners'));
    end

    requireField(cfg,'datumPointsCAD');

    if isempty(cfg.weights)
        cfg.weights = ones(size(cfg.datumPointsCAD,1),1);
    end
end


function meshInput = meshInputFromCfg(cfg,prefix)
    meshField = [prefix,'Mesh'];
    fileField = [prefix,'MeshFile'];

    if isfield(cfg,meshField) && ~isempty(cfg.(meshField))
        meshInput = cfg.(meshField);
    elseif isfield(cfg,fileField) && ~isempty(cfg.(fileField))
        meshInput = cfg.(fileField);
    else
        error("Supply cfg.%s or cfg.%s.",meshField,fileField);
    end
end


function requireField(s,name)
    if ~isfield(s,name) || isempty(s.(name))
        error("Missing required cfg field: %s",name);
    end
end


function X = normalizeRows(X)
    rowLength = vecnorm(X,2,2);
    if any(rowLength < eps)
        error("A zero-length datum normal was supplied.");
    end
    X = X./rowLength;
end


function plotDatumAlignment(TRcad,TRscan,target,measured,normals,method)
    figure('Name',['Datum alignment - ',char(method)]);
    hold on;
    axis equal;
    grid on;
    view(3);

    trisurf(TRcad.ConnectivityList, ...
        TRcad.Points(:,1),TRcad.Points(:,2),TRcad.Points(:,3), ...
        'FaceAlpha',0.12,'EdgeColor','none');

    trisurf(TRscan.ConnectivityList, ...
        TRscan.Points(:,1),TRscan.Points(:,2),TRscan.Points(:,3), ...
        'FaceAlpha',0.22,'EdgeColor','none');

    scatter3(target(:,1),target(:,2),target(:,3),60,'filled');
    scatter3(measured(:,1),measured(:,2),measured(:,3),35,'filled');

    diagonal = norm(max(TRcad.Points,[],1)-min(TRcad.Points,[],1));
    arrowScale = 0.03*diagonal;

    quiver3(target(:,1),target(:,2),target(:,3), ...
        arrowScale*normals(:,1),arrowScale*normals(:,2), ...
        arrowScale*normals(:,3),0,'LineWidth',1.2);

    for i = 1:size(target,1)
        plot3([target(i,1),measured(i,1)], ...
              [target(i,2),measured(i,2)], ...
              [target(i,3),measured(i,3)],'k-');
        text(target(i,1),target(i,2),target(i,3),sprintf('  %d',i));
    end

    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['Datum alignment: ',char(method)]);
    legend('Nominal CAD','Aligned distorted mesh', ...
        'CAD targets','Projected points','CAD normals');
end
