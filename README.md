# CT_mesh_handling

## The filtering of the meshes with meshlab:

### Before decimation:
- It might be necessary to split the mesh and clean elements separately.
  - Selection → Select non Manifold Edges
  - Selection → Delete selected vertices
  - Selection → Select non Manifold Vertices
  - Selection → Delete selected vertices
  - Quality Measure and Computation → Compute Topological Measures
  - Selection → Delete selected vertices
  - Remeshing, Simplification and Reconstruction → Close Holes (30mm used)
  - Cleaning and repairing → Remove Unreferenced Vertices
  - Quality Measure and Computation → Compute Geometric Measures
  - Quality Measure and Computation →C ompute Topological Measures


### Decimation:
- Remeshing, Simplification and Reconstruction → Simplification: Quadric edge collapse decimation
	Comments: it is only possible to choose number of faces, not vertices. A loop will be 	applied to repeat until the the given number is obtained. Can be used with Percentage reduction to get the correct number of vertices (6 relevant places!)
  - Percentage reduction: to be computed, in case_01: 0.050486
  - Quality threshold: 1.0
  - Preserve Boundary of the mesh: Off
  - Preserve Normal: On
  - Preserve Topology: On
  - Optimal position of simplified vertices: On
  - Weighted Simplification: Off
  - Planar Simplification: On
  - Post-simplification cleaning: Off
  - Simplify only selected faces: Off

###TODO:
The point correspondence must be established between meshes, then they can be applied to the atlas class.

## Remeshing:
- Clean the non manifold vertexes with selection of non manifold edges and deletion (filter doesn’t work!) then fill holes and check the geometry and topology. SAVE!
  - Use one of the subdivision filters – all except midpoint and catmull-clark (former is VTK version of it and latter builds quads)
- Results:
  - LS3 loop works better than loop
  - Setting the length to 0 yields good results
  - Enhanced regularization helps
  - It is poossible to select the faces with concrete edge lengths to run decimation on them. With a proper threshold, it could lead to nice and uniform meshes.
  - Butterfly creates patches on length 0, but it works nicely on 2
  - Use the Hausdorff distance to estimate quality and choose the best solution
  - Remeshing, Simplification and Reconstruction → Laplacian smooth (surface preserve):
  - play with angle and number of iterations

### Old (marching cubes)
- Remeshing, Simplification and Reconstruction → Uniform Mesh Resampling
  - Precision ‘world unit’: 0.5 (mm)
  - Offset ‘world unit’: 0
  - Clean Vertices: On
  - Discretize: Off
  - Multisample: On
  - Absolute Distance: Off