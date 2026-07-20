%% EXAMPLE: COMPARE RPS AND POINT-TO-POINT ALIGNMENT
clear;
clc;

load nominalSTL.mat nominalSTL
load prealignedDistortedSTL.mat distortedSTL

cfg = struct;
cfg.cadMesh = nominalSTL;
cfg.scanMesh = distortedSTL;
cfg.preAligned = true;
cfg.meshMergeTolerance = 0;

cfg.datumPointsCAD = readmatrix("wenzel_nominal_targets.csv");

% Use WENZEL vectors when available. Otherwise local CAD mesh normals are used.
% cfg.datumNormalsCAD = readmatrix("wenzel_target_vectors.csv");

cfg.maxProjectionDistance = 3.0;
cfg.comparisonSampleCount = 3000;
cfg.comparePlot = true;

comparison = compare_alignment_methods(cfg);

writetable(comparison.summary,"alignment_method_comparison.csv");
