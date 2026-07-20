%% EXAMPLE: PRE-ALIGNED STL STRUCTURES
% Your imported structures may contain:
%   nominalSTL.faces / nominalSTL.vertices
%   distortedSTL.faces / distortedSTL.vertices
%
% Repeated coordinates with different vertex IDs are accepted and merged.

clear;
clc;

% Load or construct your two STL structures here.
load nominalSTL.mat nominalSTL
load prealignedDistortedSTL.mat distortedSTL

% Optional conversion check by itself:
[TRnominal,nominalInfo] = normalize_stl_mesh(nominalSTL,0);
[TRdistorted,distortedInfo] = normalize_stl_mesh(distortedSTL,0);

fprintf("Nominal removed duplicate vertices: %d\n", ...
    nominalInfo.verticesRemoved);
fprintf("Distorted removed duplicate vertices: %d\n", ...
    distortedInfo.verticesRemoved);

cfg = struct;
cfg.cadMesh = nominalSTL;
cfg.scanMesh = distortedSTL;
cfg.preAligned = true;
cfg.meshMergeTolerance = 0;

% Coordinates exported from WENZEL:
cfg.datumPointsCAD = readmatrix("wenzel_nominal_targets.csv");

% Preferred when available:
% cfg.datumNormalsCAD = readmatrix("wenzel_target_vectors.csv");
%
% When omitted, the code derives local normals from nominalSTL.
cfg.method = "rps";
cfg.maxProjectionDistance = 3.0;
cfg.plot = true;

result = vane_rps_alignment(cfg);

% Result using the same exchange-style field names:
alignedSTL.faces = result.alignedScanMesh.faces;
alignedSTL.vertices = result.alignedScanMesh.vertices;

save alignedSTL.mat alignedSTL result
