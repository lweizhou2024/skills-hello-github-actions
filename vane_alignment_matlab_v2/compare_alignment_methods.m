function comparison = compare_alignment_methods(baseCfg)
%COMPARE_ALIGNMENT_METHODS Run RPS and point alignment from the same input.
%
% baseCfg should normally contain:
%   cadMesh
%   scanMesh                pre-aligned mesh
%   preAligned = true
%   datumPointsCAD
% Optional:
%   datumNormalsCAD
%
% The function reports datum residuals and whole-surface diagnostics for both
% methods. Do not choose a method only because it has the lowest whole-surface
% RMS: a global fit can hide twist or bow.

    if ~isfield(baseCfg,'preAligned')
        baseCfg.preAligned = true;
    end

    cfgRPS = baseCfg;
    cfgRPS.method = "rps";
    cfgRPS.plot = false;
    resultRPS = vane_rps_alignment(cfgRPS);

    cfgPoint = baseCfg;
    cfgPoint.method = "point";
    cfgPoint.plot = false;
    resultPoint = vane_rps_alignment(cfgPoint);

    if isfield(baseCfg,'comparePlot')
        comparePlot = logical(baseCfg.comparePlot);
    else
        comparePlot = true;
    end

    checkOptions = struct;
    checkOptions.plot = comparePlot;
    checkOptions.sampleCount = getOption(baseCfg,'comparisonSampleCount',2500);
    checkOptions.mergeTolerance = getOption(baseCfg,'meshMergeTolerance',0);
    checkOptions.datumPointsCAD = resultRPS.datumPointsCAD;
    checkOptions.datumNormalsCAD = resultRPS.datumNormalsCAD;
    checkOptions.maxProjectionDistance = ...
        getOption(baseCfg,'maxProjectionDistance',5);

    cadInput = getMeshInput(baseCfg,'cad');

    checkOptions.label = "RPS alignment";
    checkRPS = check_stl_alignment( ...
        cadInput,resultRPS.alignedScanMesh,checkOptions);

    checkOptions.label = "Point-to-point alignment";
    checkPoint = check_stl_alignment( ...
        cadInput,resultPoint.alignedScanMesh,checkOptions);

    Method = ["RPS";"Point-to-point"];
    AddedRotationDeg = [rotation_angle_deg(resultRPS.R); ...
                        rotation_angle_deg(resultPoint.R)];
    AddedTranslation = [norm(resultRPS.t);norm(resultPoint.t)];
    DatumNormalRMS = [resultRPS.rps.normalRMS; ...
                      resultPoint.rps.normalRMS];
    DatumPointRMS = [resultRPS.rps.pointRMS; ...
                     resultPoint.rps.pointRMS];
    SymmetricSurfaceRMS = [checkRPS.symmetricRMS; ...
                           checkPoint.symmetricRMS];
    CadToScanP95 = [checkRPS.cadToScan.p95; ...
                    checkPoint.cadToScan.p95];
    ScanToCadP95 = [checkRPS.scanToCad.p95; ...
                    checkPoint.scanToCad.p95];
    ResidualBestFitRotationDeg = [ ...
        checkRPS.residualBestFit.rotationDeg; ...
        checkPoint.residualBestFit.rotationDeg];
    ResidualBestFitTranslation = [ ...
        checkRPS.residualBestFit.translationMagnitude; ...
        checkPoint.residualBestFit.translationMagnitude];

    summary = table(Method,AddedRotationDeg,AddedTranslation, ...
        DatumNormalRMS,DatumPointRMS,SymmetricSurfaceRMS, ...
        CadToScanP95,ScanToCadP95,ResidualBestFitRotationDeg, ...
        ResidualBestFitTranslation);

    disp(summary);

    comparison.summary = summary;
    comparison.rps.result = resultRPS;
    comparison.rps.check = checkRPS;
    comparison.point.result = resultPoint;
    comparison.point.check = checkPoint;
end


function value = getOption(s,name,defaultValue)
    if isfield(s,name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = defaultValue;
    end
end


function meshInput = getMeshInput(cfg,prefix)
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
