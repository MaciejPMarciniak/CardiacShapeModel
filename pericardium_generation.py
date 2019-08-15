import vtk
import numpy as np
from vtk.util.colors import red


def get_external_surface(_mesh):
    _center = np.zeros(3)
    _bounds = np.zeros(6)
    _ray_start = np.zeros(3)
    cell_id = vtk.mutable(-1)
    xyz = np.zeros(3)
    pcoords = np.zeros(3)
    t = vtk.mutable(0)
    sub_id = vtk.mutable(0)
    _surf = 1.1

    _mesh.GetOutput().GetCenter(_center)
    _mesh.GetOutput().GetPoints().GetBounds(_bounds)
    for j in range(3):
        _ray_start[j] = _bounds[2 * j + 1] * _surf

    cell_locator = vtk.vtkCellLocator()
    cell_locator.SetDataSet(_mesh.GetOutput())
    cell_locator.BuildLocator()
    cell_locator.IntersectWithLine(_ray_start, _center, 0.0001, t, xyz, pcoords, sub_id, cell_id)

    connectivity_filter = vtk.vtkConnectivityFilter()
    connectivity_filter.SetInputConnection(_mesh.GetOutputPort())
    connectivity_filter.SetExtractionModeToCellSeededRegions()
    connectivity_filter.InitializeSeedList()
    connectivity_filter.AddSeed(cell_id)
    connectivity_filter.Update()
    return ug2pd(connectivity_filter)


def visualize(mesh):
    lineMapper = vtk.vtkPolyDataMapper()
    lineMapper.SetInputConnection(mesh.GetOutputPort())
    lineActor = vtk.vtkActor()
    lineActor.SetMapper(lineMapper)
    lineActor.GetProperty().SetColor(red)
    lineActor.GetProperty().SetOpacity(0.5)
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    ren.AddActor(lineActor)
    ren.SetBackground(1, 1, 1)
    renWin.SetSize(600, 600)
    camera = vtk.vtkCamera()
    ren.SetActiveCamera(camera)
    iren.Initialize()
    renWin.Render()
    iren.Start()


def _calculate_bounding_box_diagonal(bounds):
    return np.sqrt(np.power(bounds[0]-bounds[1], 2) +
                   np.power(bounds[2]-bounds[3], 2) +
                   np.power(bounds[4]-bounds[5], 2))


def _calculate_maximum_distance(bounds, target_offset):
    d = _calculate_bounding_box_diagonal(bounds)
    return target_offset/d


def smooth_window(mesh, number_of_iterations=15, pass_band=0.05):
    smooth = vtk.vtkWindowedSincPolyDataFilter()
    smooth.SetInputConnection(mesh.GetOutputPort())
    smooth.SetNumberOfIterations(number_of_iterations)
    smooth.BoundarySmoothingOn()
    smooth.FeatureEdgeSmoothingOff()
    smooth.SetPassBand(pass_band)
    smooth.NonManifoldSmoothingOn()
    smooth.NormalizeCoordinatesOn()
    smooth.Update()
    return smooth


def create_surrounding_surface(mesh, distance=1):
    # Create implicit model with vtkImplicitModeller.
    # The contour filter creates a surface mesh
    # at the 'distance' (in mesh's units) from the provided geometry.

    print('implicit_modeller')
    _bounds = np.zeros(6)
    bounds = np.array(mesh.GetOutput().GetPoints().GetBounds())
    max_dist = _calculate_maximum_distance(bounds, distance)
    imp = vtk.vtkImplicitModeller()
    imp.SetInputConnection(surf.GetOutputPort())
    imp.SetSampleDimensions(200, 200, 200)
    imp.SetMaximumDistance(max_dist)
    imp.SetScaleToMaximumDistance(1)
    imp.SetModelBounds(*(bounds * 1.5))
    imp.CappingOn()
    imp.SetCapValue(20.0)
    imp.Update()
    print('contouring')
    contour = vtk.vtkContourFilter()
    contour.SetInputConnection(imp.GetOutputPort())
    contour.GenerateTrianglesOn()
    contour.SetValue(0, 10.0)
    contour.Update()
    print('smoothing')
    contour = smooth_window(contour, 30)
    return get_external_surface(contour)


def append_meshes(mesh1, mesh2):
    peri = vtk.vtkAppendFilter()
    peri.AddInputConnection(mesh1.GetOutputPort())
    peri.AddInputConnection(mesh2.GetOutputPort())
    peri.Update()
    return ug2pd(peri)


def ug2pd(mesh):
    surf = vtk.vtkDataSetSurfaceFilter()
    surf.SetInputConnection(mesh.GetOutputPort())
    surf.Update()
    return surf


def save_mesh(mesh, path):
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputConnection(mesh.GetOutputPort())
    writer.SetFileName(path)
    writer.Update()
    writer.Write()
    print('{} written succesfully'.format(path))


reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("/home/mat/Deformetrica/deterministic_atlas_ct/Main_modes_export/Template_tetra.vtk")
reader.Update()

# Extract surface of the volumetric mesh
surf = vtk.vtkDataSetSurfaceFilter()
surf.SetInputConnection(reader.GetOutputPort())
surf.Update()

# Create external surface
contour = create_surrounding_surface(surf)

# Get external surface only
external_temp_surf = get_external_surface(surf)
pericardium = append_meshes(contour, external_temp_surf)

# Save results
output_filename = "/home/mat/Deformetrica/deterministic_atlas_ct/Temp/Ale.vtk"
heart_surf_filename = "/home/mat/Deformetrica/deterministic_atlas_ct/Temp/Temp_surf.vtk"
peri_surf_filename = "/home/mat/Deformetrica/deterministic_atlas_ct/Temp/peri_surf.vtk"

save_mesh(pericardium, output_filename)
save_mesh(external_temp_surf, heart_surf_filename)
save_mesh(contour, peri_surf_filename)
