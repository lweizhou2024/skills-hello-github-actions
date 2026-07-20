function report = check_stl_alignment(cadMesh,scanMesh,options)
%CHECK_STL_ALIGNMENT Quantitatively and visually compare two aligned meshes.
%
% report = check_stl_alignment(cadMesh,scanMesh,options)
%
% This function DOES NOT change the supplied meshes. It reports:
%   - CAD-to-scan closest-surface distances
%   - scan-to-CAD closest-surface distances
%   - symmetric RMS distance
%   - signed CAD-normal deviations
%   - optional RPS target residuals
%   - optional residual trimmed-ICP rigid correction
%
% options fields:
%   mergeTolerance       default 0
%   sampleCount          default 2500 points in each direction
%   neighborhoodRings    default 2
%   plot                 default true
%   label                default "Alignment check"
%   computeResidualICP   default true
%   icpIterations        default 8
%   icpKeepFraction      default 0.85
%
% Optional datum check:
%   datumPointsCAD
%   datumNormalsCAD
%   maxProjectionDistance default 5
%
% Important:
% A lower whole-surface RMS does not automatically mean a more correct datum
% alignment. A global best fit can hide physical twist or bow.

    if nargin < 3
        options = struct;
    end
    options = defaults(options);

    [TRcad,cadInfo] = normalize_stl_mesh(cadMesh,options.mergeTolerance);
    [TRscan,scanInfo] = normalize_stl_mesh(scanMesh,options.mergeTolerance);

    cadSampleID = sampleIDs(size(TRcad.Points,1),options.sampleCount);
    scanSampleID = sampleIDs(size(TRscan.Points,1),options.sampleCount);

    cadSample = TRcad.Points(cadSampleID,:);
    scanSample = TRscan.Points(scanSampleID,:);

    [scanPointForCAD,~,cadToScanDistance] = closest_points_on_mesh( ...
        cadSample,TRscan,options.neighborhoodRings);

    [cadPointForScan,~,scanToCADDistance] = closest_points_on_mesh( ...
        scanSample,TRcad,options.neighborhoodRings);

    cadVertexNormal = vertexNormal(TRcad);
    sampleNormal = cadVertexNormal(cadSampleID,:);
    sampleNormal = normalizeRowsSafe(sampleNormal);

    signedNormalDeviation = sum( ...
        (scanPointForCAD-cadSample).*sampleNormal,2);

    report.cadToScan = distanceStatistics(cadToScanDistance);
    report.scanToCad = distanceStatistics(scanToCADDistance);
    report.symmetricRMS = sqrt( ...
        (sum(cadToScanDistance.^2)+sum(scanToCADDistance.^2))/ ...
        (numel(cadToScanDistance)+numel(scanToCADDistance)));

    report.signedCadNormal = signedStatistics(signedNormalDeviation);
    report.cadSamplePoints = cadSample;
    report.scanPointsAtCadSamples = scanPointForCAD;
    report.signedNormalDeviation = signedNormalDeviation;

    report.cadMeshInfo = cadInfo;
    report.scanMeshInfo = scanInfo;

    if isfield(options,'datumPointsCAD') && ...
            ~isempty(options.datumPointsCAD)
        if ~isfield(options,'datumNormalsCAD') || ...
                isempty(options.datumNormalsCAD)
            [datumNormals,normalInfo] = cad_normals_from_mesh( ...
                options.datumPointsCAD,TRcad,1,20);
        else
            datumNormals = normalizeRowsSafe(options.datumNormalsCAD);
            normalInfo.source = "supplied";
        end

        [datumMeasured,projectionInfo] = project_targets_to_mesh( ...
            options.datumPointsCAD,datumNormals,TRscan, ...
            options.maxProjectionDistance,true);

        datumNormalResidual = sum( ...
            (datumMeasured-options.datumPointsCAD).*datumNormals,2);
        datumPointResidual = vecnorm( ...
            datumMeasured-options.datumPointsCAD,2,2);

        report.datum.pointsCAD = options.datumPointsCAD;
        report.datum.normalsCAD = datumNormals;
        report.datum.measuredPoints = datumMeasured;
        report.datum.normalResiduals = datumNormalResidual;
        report.datum.pointResiduals = datumPointResidual;
        report.datum.normalRMS = sqrt(mean(datumNormalResidual.^2));
        report.datum.pointRMS = sqrt(mean(datumPointResidual.^2));
        report.datum.normalInfo = normalInfo;
        report.datum.projectionInfo = projectionInfo;
    else
        report.datum = [];
    end

    if options.computeResidualICP
        report.residualBestFit = residualTrimmedICP( ...
            TRscan,TRcad,options.sampleCount, ...
            options.neighborhoodRings,options.icpIterations, ...
            options.icpKeepFraction);
    else
        report.residualBestFit = [];
    end

    fprintf("\n%s\n",options.label);
    fprintf("  CAD -> scan RMS        : %.6g\n",report.cadToScan.rms);
    fprintf("  CAD -> scan P95        : %.6g\n",report.cadToScan.p95);
    fprintf("  Scan -> CAD RMS        : %.6g\n",report.scanToCad.rms);
    fprintf("  Scan -> CAD P95        : %.6g\n",report.scanToCad.p95);
    fprintf("  Symmetric RMS          : %.6g\n",report.symmetricRMS);
    fprintf("  Signed-normal mean     : %.6g\n",report.signedCadNormal.mean);
    fprintf("  Signed-normal RMS      : %.6g\n",report.signedCadNormal.rms);

    if ~isempty(report.datum)
        fprintf("  Datum normal RMS       : %.6g\n",report.datum.normalRMS);
        fprintf("  Datum point RMS        : %.6g\n",report.datum.pointRMS);
    end

    if ~isempty(report.residualBestFit)
        fprintf("  Residual ICP rotation  : %.6g degrees\n", ...
            report.residualBestFit.rotationDeg);
        fprintf("  Residual ICP translation: %.6g\n", ...
            report.residualBestFit.translationMagnitude);
    end

    if options.plot
        plotOverlay(TRcad,TRscan,options.label);
        plotDeviationMap(cadSample,signedNormalDeviation,options.label);
    end
end


function options = defaults(options)
    values.mergeTolerance = 0;
    values.sampleCount = 2500;
    values.neighborhoodRings = 2;
    values.plot = true;
    values.label = "Alignment check";
    values.computeResidualICP = true;
    values.icpIterations = 8;
    values.icpKeepFraction = 0.85;
    values.maxProjectionDistance = 5;

    names = fieldnames(values);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(options,name) || isempty(options.(name))
            options.(name) = values.(name);
        end
    end
end


function IDs = sampleIDs(numberOfPoints,sampleCount)
    count = min(numberOfPoints,max(1,round(sampleCount)));
    IDs = unique(round(linspace(1,numberOfPoints,count))).';
end


function stats = distanceStatistics(values)
    values = values(:);
    stats.count = numel(values);
    stats.mean = mean(values);
    stats.median = percentile(values,50);
    stats.rms = sqrt(mean(values.^2));
    stats.p90 = percentile(values,90);
    stats.p95 = percentile(values,95);
    stats.p99 = percentile(values,99);
    stats.max = max(values);
end


function stats = signedStatistics(values)
    values = values(:);
    stats.count = numel(values);
    stats.mean = mean(values);
    stats.median = percentile(values,50);
    stats.rms = sqrt(mean(values.^2));
    stats.min = min(values);
    stats.max = max(values);
    stats.p05 = percentile(values,5);
    stats.p95 = percentile(values,95);
end


function value = percentile(values,percentage)
    sorted = sort(values(:));
    if isempty(sorted)
        value = NaN;
        return;
    end
    if numel(sorted) == 1
        value = sorted(1);
        return;
    end

    location = 1+(numel(sorted)-1)*percentage/100;
    lower = floor(location);
    upper = ceil(location);
    fraction = location-lower;
    value = sorted(lower)*(1-fraction)+sorted(upper)*fraction;
end


function output = residualTrimmedICP( ...
        movingTR,fixedTR,sampleCount,ringDepth,iterations,keepFraction)

    movingIDs = sampleIDs(size(movingTR.Points,1),sampleCount);
    originalSample = movingTR.Points(movingIDs,:);

    R = eye(3);
    t = zeros(3,1);
    historyRMS = zeros(iterations,1);

    for iteration = 1:iterations
        moved = apply_rigid(originalSample,R,t);
        [closest,~,distance] = closest_points_on_mesh( ...
            moved,fixedTR,ringDepth);

        cutoff = percentile(distance,100*keepFraction);
        keep = distance <= cutoff;

        [incrementR,incrementT,~] = weighted_kabsch( ...
            moved(keep,:),closest(keep,:),[]);

        R = incrementR*R;
        t = incrementR*t+incrementT;

        movedAfter = apply_rigid(originalSample,R,t);
        [~,~,distanceAfter] = closest_points_on_mesh( ...
            movedAfter,fixedTR,ringDepth);
        historyRMS(iteration) = sqrt(mean(distanceAfter.^2));

        if rotation_angle_deg(incrementR) < 1e-7 && norm(incrementT) < 1e-7
            historyRMS = historyRMS(1:iteration);
            break;
        end
    end

    output.R = R;
    output.t = t;
    output.rotationDeg = rotation_angle_deg(R);
    output.translationMagnitude = norm(t);
    output.historyRMS = historyRMS;
    output.finalRMS = historyRMS(end);
    output.keepFraction = keepFraction;
end


function X = normalizeRowsSafe(X)
    rowLength = vecnorm(X,2,2);
    bad = rowLength < eps;
    if any(bad)
        warning("Some vertex normals were undefined; replacing them locally.");
        X(bad,:) = repmat([0 0 1],sum(bad),1);
        rowLength = vecnorm(X,2,2);
    end
    X = X./rowLength;
end


function plotOverlay(TRcad,TRscan,label)
    figure('Name',[char(label),' - overlay']);
    hold on;
    axis equal;
    grid on;
    view(3);

    trisurf(TRcad.ConnectivityList, ...
        TRcad.Points(:,1),TRcad.Points(:,2),TRcad.Points(:,3), ...
        'FaceAlpha',0.18,'EdgeColor','none');

    trisurf(TRscan.ConnectivityList, ...
        TRscan.Points(:,1),TRscan.Points(:,2),TRscan.Points(:,3), ...
        'FaceAlpha',0.28,'EdgeColor','none');

    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title([char(label),' — mesh overlay']);
    legend('Nominal CAD','Compared mesh');
end


function plotDeviationMap(points,signedDeviation,label)
    figure('Name',[char(label),' - deviation']);
    scatter3(points(:,1),points(:,2),points(:,3), ...
        14,signedDeviation,'filled');
    axis equal;
    grid on;
    view(3);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    colorbar;
    title([char(label),' — signed CAD-normal deviation']);
end
