Vane Alignment MATLAB Package — Version 2
=========================================

What changed
------------
1. Pre-aligned mesh input is now supported directly.
2. Common File Exchange STL structs are accepted.
3. Repeated vertices are merged before creating MATLAB triangulation objects.
4. Two supplied STL meshes can be checked without changing their alignment.
5. RPS and point-to-point final alignments can be compared automatically.

Typical pre-aligned call
------------------------
cfg.cadMesh = nominalSTL;
cfg.scanMesh = distortedSTL;
cfg.preAligned = true;
cfg.datumPointsCAD = readmatrix("wenzel_nominal_targets.csv");
cfg.method = "rps";
result = vane_rps_alignment(cfg);

Supported mesh forms
--------------------
- triangulation object
- STL filename
- struct with:
    faces / vertices
    face / vertex
    ConnectivityList / Points

Duplicated STL vertices
-----------------------
A mesh such as:

vertices:
  1 = [0 0 0]
  2 = [1 0 0]
  3 = [0 1 0]
  4 = [1 0 0]   repeated coordinate with a new ID
  ...

faces:
  [1 2 3]
  [4 5 6]

is valid STL-style data. normalize_stl_mesh merges equal coordinates and updates
the face IDs. This enables reliable adjacency, vertex normals, and local surface
searches.

Use:
    [TR,info] = normalize_stl_mesh(myStruct,0);

Set a positive merge tolerance only when coordinates that should be identical
differ slightly because of rounding. Keep it much smaller than the smallest
real geometric feature.

Checking an alignment
---------------------
report = check_stl_alignment(nominalSTL,alignedSTL,options);

Main diagnostics:
- CAD-to-scan RMS and P95 distance
- scan-to-CAD RMS and P95 distance
- symmetric surface RMS
- signed CAD-normal deviation
- optional datum normal RMS
- optional residual trimmed-ICP correction

The residual ICP is diagnostic only. A large residual correction suggests that
the two meshes still contain a global rigid offset. It is not evidence that ICP
is the correct production datum alignment.

Comparing methods
-----------------
comparison = compare_alignment_methods(cfg);

The summary table compares:
- added rotation and translation
- datum normal RMS
- datum 3-D point RMS
- whole-surface symmetric RMS
- directional P95 surface deviations
- remaining global best-fit correction

Interpretation
--------------
For RPS, datum normal RMS is the primary alignment residual.
For point-to-point, datum 3-D point RMS is naturally favored.
Whole-surface RMS is useful for understanding the consequences of the
alignment, but it should not be used alone to select the datum method. A global
best fit can reduce RMS by absorbing true twist, bow, or band-to-band movement.

Normals
-------
XYZ target coordinates alone do not define normals. Preferred order:
1. Export WENZEL target vectors I,J,K or Nx,Ny,Nz.
2. Derive local normals from a fine nominal CAD mesh.
3. Manually assign directions only when they are known from the RPS definition.

Targets on sharp edges or blended junctions are ambiguous. Move them into the
interior of the intended surface or obtain the WENZEL vectors.

MATLAB compatibility
--------------------
The package uses triangulation, stlread, nearestNeighbor, vertexAttachments,
vertexNormal, fminsearch, and standard plotting functions. No Statistics
Toolbox is required.
