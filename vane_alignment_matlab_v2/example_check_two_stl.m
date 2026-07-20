%% EXAMPLE: CHECK TWO ALREADY-ALIGNED STL STRUCTURES
clear;
clc;

load nominalSTL.mat nominalSTL
load alignedCandidateSTL.mat alignedCandidateSTL

options = struct;
options.label = "Candidate alignment";
options.sampleCount = 3000;
options.plot = true;
options.computeResidualICP = true;
options.mergeTolerance = 0;

% Optional: also check WENZEL datum target residuals.
options.datumPointsCAD = readmatrix("wenzel_nominal_targets.csv");
% options.datumNormalsCAD = readmatrix("wenzel_target_vectors.csv");
options.maxProjectionDistance = 3.0;

report = check_stl_alignment( ...
    nominalSTL,alignedCandidateSTL,options);
