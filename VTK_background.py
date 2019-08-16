import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import os
import glob
from PIL import Image

# 01. LV myocardium (endo + epi)
# 02. RV myocardium (endo + epi)
# 03. LA myocardium (endo + epi)
# 04. RA myocardium (endo + epi)
#
# 05. Aorta
# 06. Pulmonary artery
#
# 07. Mitral valve
# 08. Triscupid valve
#
# 09. Aortic valve
# 10. Pulmonary valve

# 11. Appendage
# 12. Left superior pulmonary vein
# 13. Left inferior pulmonary vein
# 14. Right inferior pulmonary vein
# 15. Right superior pulmonary vein
#
# 16. Superior vena cava
# 17. Inferior vena cava

# 18. Appendage border
# 19. Right inferior pulmonary vein border
# 20. Left inferior pulmonary vein border
# 21. Left superior pulmonary vein border
# 22. Right superior pulmonary vein border
# 23. Superior vena cava border
# 24. Inferior vena cava border


class Model:

    list_of_elements = ['LV', 'RV', 'LA', 'RA', 'AO', 'PA', 'MV', 'TV', 'AV', 'PV',
                        'APP', 'LSPV', 'LIPV', 'RIPV', 'RSPV', 'SVC', 'IVC',
                        'AB', 'RIPVB', 'LIPVB', 'LSPVB', 'RSPVB', 'SVCB', 'IVCB']

    # TODO: Add the dictionary with labels, to use in alignment

    def __init__(self, filename='h_case06.vtk', to_polydata=False):

        self.filename, self.input_type = filename.split('.')
        print(self.filename)
        # Write vtk output to a file
        w = vtk.vtkFileOutputWindow()
        w.SetFileName(self.filename.split('/')[0] + '/errors.txt')
        vtk.vtkOutputWindow.SetInstance(w)

        print('Reading the data from {}.{}...'.format(self.filename, self.input_type))
        if self.input_type == 'obj':
            self.mesh, self.scalar_range = self.read_obj()
        elif self.input_type == 'vtp':
            self.mesh, self.scalar_range = self.read_vtp()
        else:
            self.mesh, self.scalar_range = self.read_vtk(to_polydata)

        # self.scalar_range = [1.0, 17.0]  # Added for particular case of CT meshes
        # print('Corrected scalar range: {}'.format(self.scalar_range))
        self.center_of_model = self.get_center(self.mesh)
        print('Model centered at: {}'.format(self.center_of_model))
        self.label = 0

    @staticmethod
    def get_center(_mesh):
        centerofmass = vtk.vtkCenterOfMass()
        centerofmass.SetInputData(_mesh.GetOutput())
        centerofmass.Update()
        return np.array(centerofmass.GetCenter())

    def visualize_mesh(self, display=True):
        # Create the mapper that corresponds the objects of the vtk file into graphics elements
        mapper = vtk.vtkDataSetMapper()
        try:
            mapper.SetInputData(self.mesh.GetOutput())
        except TypeError:
            print('Can\'t get output directly')
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(self.mesh.GetOutputPort())
        mapper.SetScalarRange(self.scalar_range)

        # Create the Actor
        camera = vtk.vtkCamera()
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        # actor.GetProperty().SetColor(vtk.util.colors.red)
        actor.GetProperty().SetOpacity(0.5)

        # Create the Renderer
        renderer = vtk.vtkRenderer()
        renderer.ResetCameraClippingRange()
        renderer.AddActor(actor)  # More actors can be added
        renderer.SetActiveCamera(camera)
        renderer.SetBackground(1, 1, 1)  # Set background to white

        # Create the RendererWindow
        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer(renderer)
        render_window.SetSize(600, 600)
        render_window.Render()

        # Display the mesh
        # noinspection PyArgumentList
        if display:
            interactor = vtk.vtkRenderWindowInteractor()
            interactor.SetRenderWindow(render_window)
            interactor.Initialize()
            interactor.Start()
        else:
            return render_window

    # -----3D rigid transformations---------------------------------------------------------------------------

    def rotate(self, alpha=0, beta=0, gamma=0, rotation_matrix=None):
        print('rotating')
        rotate = vtk.vtkTransform()
        if rotation_matrix is not None:
            translation_matrix = np.eye(4)
            translation_matrix[:-1, :-1] = rotation_matrix
            print('Translation matrix (rotation):\n', translation_matrix)
            rotate.SetMatrix(translation_matrix.ravel())
        else:
            rotate.Identity()
            rotate.RotateX(alpha)
            rotate.RotateY(beta)
            rotate.RotateZ(gamma)
        transformer = vtk.vtkTransformFilter()
        transformer.SetInputConnection(self.mesh.GetOutputPort())
        transformer.SetTransform(rotate)
        transformer.Update()
        self.mesh = transformer
        self.center_of_model = self.get_center(self.mesh)

    def scale(self, factor=(0.001, 0.001, 0.001)):
        print('scaling')
        scale = vtk.vtkTransform()
        scale.Scale(factor[0], factor[1], factor[2])
        transformer = vtk.vtkTransformFilter()
        transformer.SetInputConnection(self.mesh.GetOutputPort())
        transformer.SetTransform(scale)
        transformer.Update()
        self.mesh = transformer
        self.center_of_model = self.get_center(self.mesh)
        print(self.center_of_model)

    def translate(self, rotation_matrix, translation_vector):
        print('translating')
        translate = vtk.vtkTransform()
        translation_matrix = np.eye(4)
        translation_matrix[:-1, :-1] = rotation_matrix
        translation_matrix[:-1, -1] = translation_vector
        print('Translation matrix:\n', translation_matrix)
        translate.SetMatrix(translation_matrix.ravel())
        transformer = vtk.vtkTransformFilter()
        transformer.SetInputConnection(self.mesh.GetOutputPort())
        transformer.SetTransform(translate)
        transformer.Update()
        self.mesh = transformer
        self.center_of_model = self.get_center(self.mesh)

    def translate_to_center(self, label=None):
        # vtkTransform.SetMatrix - enables for applying 4x4 transformation matrix to the meshes
        # if label is provided, translates to the center of the element with that label
        print('translating o center')
        translate = vtk.vtkTransform()
        if label is not None:
            central_element = self.threshold(label, label)
            center_of_element = self.get_center(central_element)
            translate.Translate(-center_of_element[0], -center_of_element[1], -center_of_element[2])
        else:
            translate.Translate(-self.center_of_model[0], -self.center_of_model[1], -self.center_of_model[2])
        translate.Update()
        transformer = vtk.vtkTransformFilter()
        transformer.SetInputConnection(self.mesh.GetOutputPort())
        transformer.SetTransform(translate)
        transformer.Update()
        self.mesh = transformer
        self.center_of_model = self.get_center(self.mesh)
        print(self.center_of_model)

    # -----Mesh manipulation----------------------------------------------------------------------------------
    def align_slice(self, a, b, c):
        print('aligning slice')
        center = np.mean((a, b, c), axis=0)
        self.translate(rotation_matrix=np.eye(3), translation_vector=-center)

        a2, b2, c2 = [x - center for x in [a, b, c]]
        _normal = calculate_plane_normal(a2, b2, c2)
        rot1 = calculate_rotation(np.array([0, 0, 1]), _normal)
        a3, b3, c3 = [rot1 @ x for x in [a2, b2, c2]]
        rot2 = calculate_rotation(np.array([0, 1, 0]), b3/np.linalg.norm(b3))
        rot3 = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        rot = rot3 @ rot2 @ rot1
        self.rotate(rotation_matrix=rot)

    def apply_modes(self, modes_with_scales):
        print('applying modes')
        for mode, scale in modes_with_scales.items():
            print('Applying ' + mode + ' multiplied by ' + str(scale))
            self.mesh.GetOutput().GetPointData().SetActiveVectors(mode)
            warp_vector = vtk.vtkWarpVector()
            warp_vector.SetInputConnection(self.mesh.GetOutputPort())
            warp_vector.SetScaleFactor(scale)
            warp_vector.Update()
            self.mesh = warp_vector

    def build_tag(self, label):
        print('building tag')
        self.label = label
        tag = vtk.vtkIdFilter()
        tag.CellIdsOn()
        tag.PointIdsOff()
        tag.SetInputConnection(self.mesh.GetOutputPort())
        tag.SetIdsArrayName('elemTag')
        tag.Update()
        self.mesh = tag

    @staticmethod
    def calculate_bounding_box_diagonal(bounds):
        return np.sqrt(np.power(bounds[0] - bounds[1], 2) +
                       np.power(bounds[2] - bounds[3], 2) +
                       np.power(bounds[4] - bounds[5], 2))

    def calculate_maximum_distance(self, bounds, target_offset):
        d = self.calculate_bounding_box_diagonal(bounds)
        return target_offset / d

    def change_tag_label(self):
        print('changing tag label')
        size = self.mesh.GetOutput().GetAttributes(1).GetArray(0).GetSize()
        for id in range(size):
            self.mesh.GetOutput().GetAttributes(1).GetArray(0).SetTuple(id, (float(self.label),))

    def clean_polydata(self, tolerance=0.005, remove_lines=False):
        print('cleaning polydata')
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputConnection(self.mesh.GetOutputPort())
        cleaner.SetTolerance(tolerance)
        cleaner.ConvertLinesToPointsOn()
        cleaner.ConvertPolysToLinesOn()
        cleaner.ConvertStripsToPolysOn()
        cleaner.Update()
        self.mesh = cleaner
        if remove_lines:
            self.mesh.GetOutput().SetLines(vtk.vtkCellArray())

    def contouring(self):
        print('contouring')
        contour = vtk.vtkContourFilter()
        contour.SetInputConnection(self.mesh.GetOutputPort())
        contour.GenerateTrianglesOn()
        contour.SetValue(0, 10.0)
        contour.Update()
        self.mesh = contour

    def decimation(self, reduction=50):
        print('decimating')
        decimation = vtk.vtkQuadricDecimation()
        decimation.SetInputConnection(self.mesh.GetOutputPort())
        decimation.VolumePreservationOn()
        decimation.SetTargetReduction(reduction / 100)  # percent of kept triangles
        decimation.Update()
        self.mesh = decimation

    def delaunay2d(self):
        print('triangulating 2D')
        delaunay2d = vtk.vtkDelaunay2D()
        delaunay2d.SetInputConnection(self.mesh.GetOutputPort())
        delaunay2d.Update()
        self.mesh = delaunay2d

    def delaunay3d(self):
        print('triangulating 3D')
        delaunay3d = vtk.vtkDelaunay3D()
        delaunay3d.SetInputConnection(self.mesh.GetOutputPort())
        delaunay3d.Update()
        self.mesh = delaunay3d

    def extract_surface(self):
        print('extracting surface')
        # Get surface of the mesh
        surface_filter = vtk.vtkDataSetSurfaceFilter()
        surface_filter.SetInputData(self.mesh.GetOutput())
        surface_filter.Update()
        self.mesh = surface_filter

    def fill_holes(self, hole_size=10.0):
        print('filling holes')
        filling_filter = vtk.vtkFillHolesFilter()
        filling_filter.SetInputConnection(self.mesh.GetOutputPort())
        filling_filter.SetHoleSize(hole_size)
        filling_filter.Update()
        self.mesh = filling_filter

    def get_external_surface(self):
        print('getting external surface')
        _center = np.zeros(3)
        _bounds = np.zeros(6)
        _ray_start = np.zeros(3)
        cell_id = vtk.mutable(-1)
        xyz = np.zeros(3)
        pcoords = np.zeros(3)
        t = vtk.mutable(0)
        sub_id = vtk.mutable(0)
        _surf = 1.1

        self.mesh.GetOutput().GetCenter(_center)
        self.mesh.GetOutput().GetPoints().GetBounds(_bounds)
        for j in range(3):
            _ray_start[j] = _bounds[2 * j + 1] * _surf

        cell_locator = vtk.vtkCellLocator()
        cell_locator.SetDataSet(self.mesh.GetOutput())
        cell_locator.BuildLocator()
        cell_locator.IntersectWithLine(_ray_start, _center, 0.0001, t, xyz, pcoords, sub_id, cell_id)

        connectivity_filter = vtk.vtkConnectivityFilter()
        connectivity_filter.SetInputConnection(self.mesh.GetOutputPort())
        connectivity_filter.SetExtractionModeToCellSeededRegions()
        connectivity_filter.InitializeSeedList()
        connectivity_filter.AddSeed(cell_id)
        connectivity_filter.Update()
        self.mesh = connectivity_filter  # UnstructuredGrid

    def implicit_modeller(self, distance):
        print('implicit modelling')
        # Create implicit model with vtkImplicitModeller at the 'distance' (in mesh's units) from the provided geometry.
        bounds = np.array(self.mesh.GetOutput().GetPoints().GetBounds())
        max_dist = self.calculate_maximum_distance(bounds, distance)
        imp = vtk.vtkImplicitModeller()
        imp.SetInputConnection(self.mesh.GetOutputPort())
        imp.SetSampleDimensions(400, 400, 400)
        imp.SetMaximumDistance(max_dist)
        imp.ScaleToMaximumDistanceOn()
        imp.SetModelBounds(*(bounds * 1.5))
        imp.CappingOn()
        imp.SetCapValue(255)
        imp.Update()
        self.mesh = imp

    def measure_average_edge_length(self):
        size = vtk.vtkCellSizeFilter()
        size.SetInputConnection(self.mesh.GetOutputPort())
        size.Update()
        print(size)

    def normals(self):
        print('getting normals')
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(self.mesh.GetOutputPort())
        normals.FlipNormalsOn()
        normals.Update()
        self.mesh = normals

    def pass_array(self):
        print('passing arrays')
        passer = vtk.vtkPassArrays()
        passer.SetInputConnection(self.mesh.GetOutputPort())
        passer.AddCellDataArray('elemTag')
        passer.Update()
        self.mesh = passer

    def resample_to_image(self, label_name='elemTag'):
        print('resampling to image')
        resampler = vtk.vtkResampleToImage()
        resampler.SetInputConnection(self.mesh.GetOutputPort())
        resampler.UseInputBoundsOff()
        bounds = np.array(self.mesh.GetOutput().GetBounds())
        bounds[:4] = bounds[:4] + 0.1 * bounds[:4]
        assert np.sum(bounds[4:] < 0.001), 'The provided slice must be 2D and must be projected on the XY plane'

        resampler.SetSamplingBounds(*bounds[:5], 1.01)
        resampler.SetSamplingDimensions(1024, 1024, 1)
        resampler.Update()

        img_as_array = vtk_to_numpy(resampler.GetOutput().GetPointData().GetArray(label_name))
        img_as_array = img_as_array.reshape((int(np.sqrt(img_as_array.shape[0])), int(np.sqrt(img_as_array.shape[0]))))

        return img_as_array

    def slice_extraction(self, origin, normal):
        print('extracting slices')
        # create a plane to cut (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
        plane = vtk.vtkPlane()
        plane.SetOrigin(*origin)
        plane.SetNormal(*normal)

        # create cutter
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(plane)
        cutter.SetInputConnection(self.mesh.GetOutputPort())
        cutter.Update()

        self.mesh = cutter

    def smooth_laplacian(self, number_of_iterations=50):
        print('laplacian smoothing')
        smooth = vtk.vtkSmoothPolyDataFilter()
        smooth.SetInputConnection(self.mesh.GetOutputPort())
        smooth.SetNumberOfIterations(number_of_iterations)
        smooth.FeatureEdgeSmoothingOff()
        smooth.BoundarySmoothingOn()
        smooth.Update()
        self.mesh = smooth

    def smooth_window(self, number_of_iterations=30, pass_band=0.05):
        print('window smoothing')
        smooth = vtk.vtkWindowedSincPolyDataFilter()
        smooth.SetInputConnection(self.mesh.GetOutputPort())
        smooth.SetNumberOfIterations(number_of_iterations)
        smooth.BoundarySmoothingOn()
        smooth.FeatureEdgeSmoothingOff()
        smooth.SetPassBand(pass_band)
        smooth.NonManifoldSmoothingOn()
        smooth.NormalizeCoordinatesOn()
        smooth.Update()
        self.mesh = smooth

    def subdivision(self, number_of_subdivisions=3):
        print('subdividing')
        self.normals()
        subdivision = vtk.vtkLinearSubdivisionFilter()
        subdivision.SetNumberOfSubdivisions(number_of_subdivisions)
        subdivision.SetInputConnection(self.mesh.GetOutputPort())
        subdivision.Update()
        self.mesh = subdivision
        self.visualize_mesh(True)

    def tetrahedralize(self, leave_tetra_only=True):
        print('creating tetrahedrons')
        tetra = vtk.vtkDataSetTriangleFilter()
        if leave_tetra_only:
            tetra.TetrahedraOnlyOn()
        tetra.SetInputConnection(self.mesh.GetOutputPort())
        tetra.Update()
        self.mesh = tetra

    def threshold(self, low=0, high=100):
        print('thresholding')
        threshold = vtk.vtkThreshold()
        threshold.SetInputConnection(self.mesh.GetOutputPort())
        threshold.ThresholdBetween(low, high)
        threshold.Update()
        # choose scalars???
        return threshold

    def ug_geometry(self):
        print('setting unstructured grid geometry')
        geometry = vtk.vtkUnstructuredGridGeometryFilter()
        print(geometry.GetDuplicateGhostCellClipping())
        geometry.SetInputConnection(self.mesh.GetOutputPort())
        geometry.Update()
        self.mesh = geometry

    def unstructured_grid_to_poly_data(self):
        print('transforming UG into PD')
        surface = vtk.vtkDataSetSurfaceFilter()
        surface.SetInputConnection(self.mesh.GetOutputPort())
        surface.Update()
        return surface

    # -----MeshInformation------------------------------------------------------------------------------------
    def get_volume(self):
        mass = vtk.vtkMassProperties()
        mass.SetInputConnection(self.mesh.GetOutputPort())
        return mass.GetVolume()

    def print_numbers(self):
        _mesh = self.mesh.GetOutput()
        print('Number of vertices: {}'.format(_mesh.GetNumberOfVerts()))
        print('Number of lines: {}'.format(_mesh.GetNumberOfLines()))
        print('Number of strips: {}'.format(_mesh.GetNumberOfStrips()))
        print('Number of polys: {}'.format(_mesh.GetNumberOfPolys()))
        print('Number of cells: {}'.format(_mesh.GetNumberOfCells()))
        print('Number of points: {}'.format(_mesh.GetNumberOfPoints()))

    # -----InputOutput----------------------------------------------------------------------------------------

    # -----Readers

    def read_vtk(self, to_polydata=False):
        # Read the source file.
        assert os.path.isfile('.' .join([self.filename, self.input_type])), \
            'File {} does not exist!'.format('.' .join([self.filename, self.input_type]))
        reader = vtk.vtkDataReader()
        reader.SetFileName('.' .join([self.filename, self.input_type]))
        reader.Update()
        print('Case ID : {}, input type: {}'.format(self.filename, self.input_type))
        if reader.IsFileUnstructuredGrid():
            print('Reading Unstructured Grid...')
            reader = vtk.vtkUnstructuredGridReader()
        elif reader.IsFilePolyData():
            print('Reading Polygonal Mesh...')
            reader = vtk.vtkPolyDataReader()
        elif reader.IsFileStructuredGrid():
            print('Reading Structured Grid...')
            reader = vtk.vtkStructuredGridReader()
        elif reader.IsFileStructuredPoints():
            print('Reading Structured Points...')
            reader = vtk.vtkStructuredPointsReader()
        elif reader.IsFileRectilinearGrid():
            print('Reading Rectilinear Grid...')
            reader = vtk.vtkRectilinearGridReader()
        else:
            print('Data format unknown...')
        reader.SetFileName(self.filename + '.' + self.input_type)
        reader.Update()  # Needed because of GetScalarRange
        scalar_range = reader.GetOutput().GetScalarRange()
        if to_polydata and not reader.IsFilePolyData():
            print('Transform to Polygonal Mesh')
            reader = self.unstructured_grid_to_poly_data(reader)
        print('Scalar range: \n{}'.format(scalar_range))
        return reader, scalar_range

    def read_vtp(self):
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName('.' .join([self.filename, self.input_type]))
        reader.Update()
        scalar_range = reader.GetOutput().GetScalarRange()
        return reader, scalar_range

    def read_obj(self):
        reader = vtk.vtkOBJReader()
        reader.SetFileName('.' .join([self.filename, self.input_type]))
        reader.Update()
        scalar_range = reader.GetOutput().GetScalarRange()
        return reader, scalar_range

    # -----Writers

    def write_mha(self):

        output_filename = self.filename + '.mha'
        # output_filename_raw = self.filename + '.raw'
        print('writing mha')

        mha_writer = vtk.vtkMetaImageWriter()
        mha_writer.SetInputConnection(self.mesh.GetOutputPort())
        mha_writer.SetFileName(output_filename)
        # mha_writer.SetRAWFileName(output_filename_raw)
        mha_writer.Write()

    def write_stl(self):
        output_filename = self.filename + '.stl'

        # Get surface of the mesh
        print('Extracting surface to save as .STL file...')
        # self.extract_surface()

        # Write file to .stl format
        stl_writer = vtk.vtkSTLWriter()
        stl_writer.SetFileName(output_filename)
        stl_writer.SetInputConnection(self.mesh.GetOutputPort())
        stl_writer.Write()
        print('{} written succesfully'.format(output_filename))

    def write_obj(self, postscript=''):
        output_filename = self.filename
        render_window = self.visualize_mesh(False)

        print('Saving PolyData in the OBJ file...')
        obj_writer = vtk.vtkOBJExporter()
        obj_writer.SetRenderWindow(render_window)
        obj_writer.SetFilePrefix(output_filename + postscript)
        obj_writer.Write()
        print('{} written succesfully'.format(output_filename + postscript + '.obj'))

    def write_png(self, postscript=''):

        print('Saving slice in PNG file...')
        output_filename = self.filename + postscript + '.png'
        image = Image.fromarray(self.resample_to_image())
        image = image.convert('L')
        image.save(output_filename, 'PNG')
        print('{} written succesfully'.format(output_filename))

    def write_vtk(self, postscript='_new', type_='PolyData'):
        output_filename = self.filename + postscript + '.vtk'
        writer = None
        if type_ == 'PolyData':
            print('Saving PolyData...')
            self.extract_surface()
            writer = vtk.vtkPolyDataWriter()
        elif type_ == 'UG':
            print('Saving Unstructured Grid...')
            writer = vtk.vtkUnstructuredGridWriter()
        else:
            exit("Select \'Polydata\' or \'UG\' as type of the saved mesh")
        writer.SetInputConnection(self.mesh.GetOutputPort())
        writer.SetFileName(output_filename)
        writer.Update()
        writer.Write()
        print('{} written succesfully'.format(output_filename))

    def write_vtk_points(self, postscript='_points'):
        output_filename = self.filename + postscript + '.vtk'

        point_cloud = vtk.vtkPolyData()
        point_cloud.SetPoints(self.mesh.GetOutput().GetPoints())
        writer = vtk.vtkPolyDataWriter()
        writer.SetInputData(point_cloud)
        writer.SetFileName(output_filename)
        writer.Update()
        writer.Write()
    # --------------------------------------------------------------------------------------------------------


# -----ApplyToChambers----------------------------------------------------------------------------------------
def apply_to_chambers(input_base, version, chambers_available, start=1, end=20, ext='_new', function_=None):
    for chamber in chambers_available:
        version += chamber
        apply_function_to_all(input_base, version, start, end, ext, function_)
# ------------------------------------------------------------------------------------------------------------


# -----Solve decimation---------------------------------------------------------------------------------------
def hill_climb_with_magnitude(initial_x, eval_function, goal, this_model, step, max_iter=500):
    best_x = initial_x
    iteration = 0
    current_result = eval_function(best_x, this_model)
    print('goal: {}'.format(goal))
    print('initial_result: {}'.format(current_result))
    while current_result != goal and iteration < max_iter:
        _step = (current_result - goal) * step
        print('step: {}'.format(_step))
        best_x += _step
        current_result = eval_function(best_x, this_model)
        iteration += 1
        print('Iteration {}, decimated [%]: {} and point count: {}:'.format(iteration, best_x, current_result))
    return best_x


def eval_decimation(_reduction, _model):
    deci = vtk.vtkQuadricDecimation()
    deci.SetInputConnection(_model.mesh.GetOutputPort())
    deci.VolumePreservationOn()
    deci.SetTargetReduction(_reduction / 100)  # percent of removed triangles
    deci.Update()
    return deci.GetOutput().GetNumberOfPoints()


def decimate_heart(full_model, case=None):
    all_points = full_model.mesh.GetOutput().GetNumberOfPoints()
    final_points = 10000
    reduction = (1 - final_points / all_points) * 100
    print('Decimating full heart...')
    _model = full_model
    _model.smooth_laplacian(200)
    reduction = hill_climb_with_magnitude(initial_x=reduction, eval_function=eval_decimation,
                                          goal=final_points, step=0.0005, this_model=_model)
    _model.decimation(reduction=reduction)
    _model.write_vtk(postscript='_decimated')
    return full_model


def decimate_chambers(full_model, case=None):
    all_points = full_model.mesh.GetOutput().GetNumberOfPoints()
    final_points = 10000
    built_points = 0
    reduction = (1 - final_points / all_points) * 100
    chambers = ['lv', 'rv', 'ra', 'la']
    case_split = case.split('.')
    sum_chamber_points = 0
    for chamber in chambers:
        print('decimating ' + chamber)
        _model = Model(case_split[0].strip('ful') + chamber + '.' + case_split[1])
        _model.smooth_laplacian(500)
        chamber_points = _model.mesh.GetOutput().GetNumberOfPoints()
        expected_number_of_points = np.round((1 - reduction / 100) * chamber_points)
        built_points += expected_number_of_points
        if chamber == 'la':
            expected_number_of_points += final_points - built_points
        reduction = hill_climb_with_magnitude(initial_x=reduction, eval_function=eval_decimation,
                                              goal=expected_number_of_points, step=0.001, this_model=_model)
        _model.decimation(reduction=reduction)
        sum_chamber_points += _model.mesh.GetOutput().GetNumberOfPoints()
        _model.normals()
        _model.write_vtk(postscript='_decimated')
        # exit('checking decimation in decimate_chambers!')
    if sum_chamber_points != final_points:
        exit('Chamber_points: {}, but should be {}'.format(sum_chamber_points, final_points))

    return full_model
# ------------------------------------------------------------------------------------------------------------


# -----Splitting----------------------------------------------------------------------------------------------
def split_chambers(_model, return_as_surface=False, return_elements=True):
    # _model.translate_to_center()

    surfaces = []
    for i in range(1, int(_model.scalar_range[1]) + 1):

        x = _model.threshold(i, i)
        surfaces.append(x)

    full_model_appended = vtk.vtkAppendFilter()
    _model.filename = os.path.join(_model.filename)
    for surf, elem in zip(surfaces, _model.list_of_elements):
        print(elem)
        if return_elements:
            _model.mesh = surf
            _model.extract_surface()
            _model.write_vtk(postscript='_'+elem)

        full_model_appended.AddInputConnection(surf.GetOutputPort())

    full_model_appended.Update()
    _model.mesh = full_model_appended
    if return_as_surface:
        # _model.translate_to_center()
        _model.extract_surface()
        _model.write_vtk(postscript='surf')
    else:
        _model.write_vtk(postscript='tetra')
    return _model


# def split_and_combine_chambers(_model, case=None):
#     # Thresholds related to arteries commented, as well as saving particular meshes
#     _model.translate_to_center()
#     _mesh = _model.mesh
#     # labels of chambers
#     lv = [1, 7, 9]
#     rv = [2, 8, 10]
#     la = [3, 11, 12, 13, 14, 15]
#     ra = [4, 16, 17]
#     aorta = [5]
#     pa = [6]
#
#     surfaces = {'surface_lv': [], 'surface_rv': [], 'surface_la': [], 'surface_ra': [],
#                 'surface_aorta': [], 'surface_pa': []}
#
#     for i in range(1, int(_model.scalar_range[1])+1):
#         surface = _model.threshold(i, i)
#         if i in lv:
#             surfaces['surface_lv'].append(surface)
#         elif i in rv:
#             surfaces['surface_rv'].append(surface)
#         elif i in la:
#             surfaces['surface_la'].append(surface)
#         elif i in ra:
#             surfaces['surface_ra'].append(surface)
#         elif i in aorta:
#             surfaces['surface_aorta'].append(surface)
#         elif i in pa:
#             surfaces['surface_pa'].append(surface)
#         else:
#             print('Label {} doesn\'t belong to any chamber of interest'.format(i))
#
#     full_model_appended = vtk.vtkAppendPolyData()
#     for surf, part in surfaces.items():
#         app = vtk.vtkAppendFilter()
#         for element in part:
#             app.AddInputConnection(element.GetOutputPort())
#         app.Update()
#         _model.mesh = app
#         _model.extract_surface()
#         full_model_appended.AddInputConnection(_model.mesh.GetOutputPort())
#         # _model.visualize_mesh()
#         # _model.write_vtk(postscript=surf)
#     full_model_appended.Update()
#     _model.mesh = full_model_appended
#
#     return _model


def change_downloaded_files_names(path='h_case_', key='surfmesh', ext='vtk'):
    files = glob.glob(os.path.join(path, '*'+key+'*'+ext))
    for i, old_file in enumerate(files):
        new_file = old_file.split('.')[0]
        print(new_file)
        os.rename(old_file, new_file+'.'+ext)
    print(files)


def change_elem_tag(_mesh, label):
    size = _mesh.GetOutput().GetAttributes(1).GetArray(0).GetSize()
    for i in range(size):
        _mesh.GetOutput().GetAttributes(1).GetArray(0).SetTuple(i, (float(label),))
    return _mesh


def assign_tags(_mesh, label_and_range_tuple=({},)):
    _mesh.GetOutput().GetAttributes(1).GetArray(0).SetName('elemTag')
    _mesh.GetOutput().GetAttributes(0).RemoveArray('elemTag')  # remove point attribute
    for label_and_range in label_and_range_tuple:
        label = label_and_range['label']
        range_of_points = label_and_range['range']
        print('Assiging label {} to {} points'.format(label, range_of_points[1]))
        for id in range(*range_of_points):
            _mesh.GetOutput().GetAttributes(1).GetArray('elemTag').SetTuple(id, (float(label),))
    return _mesh


def remove_array(_model, array_name):

    _model.mesh.GetOutput().GetAttributes(1).RemoveArray(array_name)  # remove point attribute
    return _model


def merge_elements(elem1, elem2):
    """
    Appends elements and returns the single connected mesh. The points in the same position in 3D are merged into one.
    :param elem1: Single element. The order of the elements pays no role.
    :param elem2: Single element.
    :return: Merged element as filter.
    """
    merger = vtk.vtkAppendFilter()
    merger.MergePointsOn()
    merger.AddInputConnection(elem1.GetOutputPort())
    merger.AddInputConnection(elem2.GetOutputPort())
    merger.Update()
    return merger


def calculate_rotation(reference_vector, target_vector):
    """
    Calculates the rotation matrix which rotates the object to align the target vector direction to reference
    vector direction. Assumes that both vectors are anchored at the beginning of the coordinate system
    :param reference_vector: Vector with referential direction. The rotation matrix will align the target_vector's
    direction to this one.
    :param target_vector:  Vector pointing to a  structure corresponding to the referential vector.
    :return: 3x3 rotation matrix (rot), where [rot @ target_vector = reference_vector] in terms of direction.
    """

    unit_reference_vector = reference_vector / np.linalg.norm(reference_vector)
    unit_target_vector = target_vector / np.linalg.norm(target_vector)
    c = unit_target_vector @ unit_reference_vector
    if c == 1:
        return np.eye(3)
    elif c == -1:
        return -np.eye(3)
    else:
        v = np.cross(unit_target_vector, unit_reference_vector)
        vx = np.array(([0,     -v[2],   v[1]],
                       [v[2],   0,     -v[0]],
                       [-v[1],  v[0],   0]))
        vx2 = vx @ vx
        return np.eye(3) + vx + vx2 / (1 + c)


def calculate_plane_normal(a, b, c):
    """
    :param a: 3D point
    :param b: 3D point
    :param c: 3D point
    :return: Vector normal to a plane which crosses the abc points.
    """
    x = np.cross(b-a, b-c)
    return x/np.linalg.norm(x)


def get_centers(_model, _labels):

    centers = []
    for lab in _labels:
        centers.append(_model.get_center(_model.threshold(lab, lab)))
    print('centers: {}'.format(centers))
    return centers


def get_translation_vector(target_markers, reference_markers):

    target_center = np.mean([target_markers[0], target_markers[1]], axis=0)
    reference_center = np.mean([reference_markers[0], reference_markers[1]], axis=0)
    return reference_center - target_center


def get_plane_alignment_rotation_matrix(target_markers, reference_markers):

    target_center = np.mean([target_markers[0], target_markers[1]], axis=0)
    reference_center = np.mean([reference_markers[0], reference_markers[1]], axis=0)
    assert np.all(target_center - reference_center < 0.001), 'The models are not position-aligned'

    target_plane_normal = calculate_plane_normal(target_markers[2], target_center, target_markers[1])
    reference_plane_normal = calculate_plane_normal(reference_markers[2], reference_center, reference_markers[1])
    return calculate_rotation(reference_plane_normal, target_plane_normal)


def get_lowest_septal_point(_model):

    lv = _model.threshold(1, 1)
    rv = _model.threshold(2, 2)
    print('---------------------------------------------------------------------')
    lv_points = set(tuple(map(tuple, vtk_to_numpy(lv.GetOutput().GetPoints().GetData()))))
    rv_points = set(tuple(map(tuple, vtk_to_numpy(rv.GetOutput().GetPoints().GetData()))))
    common_points = lv_points.intersection(rv_points)

    valve_centers = get_centers(_model, (7, 8))
    center = np.mean([valve_centers[0], valve_centers[1]], axis=0)
    norms = np.array([[x, np.linalg.norm((center - x))] for x in common_points])
    lowest_septal_point = np.array(norms[norms[:, 1] == np.max(norms[:, 1])][0, 0])
    return lowest_septal_point


def get_vector_alignment_rotation_matrix(target_markers, reference_markers):

    target_center = np.mean([target_markers[0], target_markers[1]], axis=0)
    reference_center = np.mean([reference_markers[0], reference_markers[1]], axis=0)
    return calculate_rotation(reference_markers[2] - reference_center, target_markers[2] - target_center)


def alignment(target_model, reference_model, labels=(7, 8)):

    if len(labels) != 2:
        exit('Check the provided labels, only 2 allowed')
    else:
        tar_markers = get_centers(target_model, labels)
        ref_markers = get_centers(reference_model, labels)
        tranlsation = get_translation_vector(tar_markers, ref_markers)  # requires 2 markers, finds middle between them
        target_model.translate(rotation_matrix=np.eye(3), translation_vector=tranlsation)
        tar_lsp = get_lowest_septal_point(target_model)
        ref_lsp = get_lowest_septal_point(reference_model)
        new_tar_markers = get_centers(target_model, labels)
        new_tar_markers.append(tar_lsp)
        ref_markers.append(ref_lsp)
        rotation1 = get_plane_alignment_rotation_matrix(new_tar_markers, ref_markers)
        rotation2 = get_vector_alignment_rotation_matrix(new_tar_markers, ref_markers)
        rotation = rotation2 @ rotation1
        target_model.translate(rotation_matrix=rotation, translation_vector=np.zeros(3))
        return target_model


def get_plax_landmarks(_model):

    lv = _model.threshold(1, 1)
    lv_points = vtk_to_numpy(lv.GetOutput().GetPoints().GetData())
    valve_centers = get_centers(_model, (7, 9))
    apex_id = np.argmax([np.linalg.norm(valve_centers[0] - x) for x in lv_points])
    apex = lv_points[apex_id]

    normal = calculate_plane_normal(*valve_centers, apex)
    origin = np.mean(np.array((apex, *valve_centers)), axis=0)

    return origin, normal, (*valve_centers, apex)


def create_plax_slices(_model):
    origin, normal, landmarks = get_plax_landmarks(_model)
    _model.slice_extraction(origin, normal)
    _model.align_slice(landmarks[2], landmarks[1], landmarks[0])
    _model.rotate(gamma=-90)
    return _model

# TODO: Make all of the functions and parameters below into a class!!!
# -----ApplyToCohort------------------------------------------------------------------------------------------


def apply_single_transformation_to_all(path, input_base, version, start=0, end=0, ext='_new', ext_type='PolyData',
                                       function_=None, args='()'):
    if function_ is not None:
        if start == end:
            cases = [os.path.join(path, f) for f in os.listdir(path) if f[-4:] == ".vtk"]
        else:
            cases = [path + '/' + input_base + str(case_no).zfill(2) + version + '.vtk' for case_no in
                     range(start, end + 1)]
        print('Cases: {}'.format(cases))
        for case in cases:
            single_model = Model(case)
            print('Executing single_model.' + function_ + args)
            exec('single_model.' + function_ + args)
            if ext is not None:
                single_model.write_vtk(postscript=ext, type_=ext_type)


def apply_function_to_all(path, input_base, version, start=1, end=20, ext='_new', ext_type='PolyData',
                          function_=None, args=''):
    if function_ is not None:
        if start == end:
            cases = [os.path.join(path, f) for f in os.listdir(path) if f[-4:] == ".vtk"]
        else:
            cases = [path + '/' + input_base + str(case_no).zfill(2) + version + '.vtk' for case_no in
                     range(start, end + 1)]
        for c, case in enumerate(cases):
            print(c, case)
            single_model = Model(case)
            if function_ == 'align_with_rotation_only' and c == 0:
                anchoring_element = single_model.get_center(single_model.threshold(7, 7))
                direction_vector = single_model.get_center(single_model.threshold(8, 8)) - anchoring_element
                target_plane_norm = single_model.get_center(single_model.threshold(9, 9))  # 9 ~ Aortic valve
                plane_norm = calculate_plane_normal(anchoring_element, direction_vector/2, target_plane_norm)
                args = args+'direction_vector=direction_vector, plane_norm=plane_norm'
            if function_ == 'alignment' and c == 0:
                reference_model = single_model
                args = args+'reference_model = reference_model'
            exec('single_model = ' + function_ + '(single_model,' + args + ')')
            if ext is not None:
                single_model.write_vtk(ext, type_=ext_type)
# ------------------------------------------------------------------------------------------------------------


def h_case_pipeline(start_=1, end_=19, path=None):
    # TODO: Problem with this solution is that the files are read and written at every step. Should become
    # TODO: a pipeline, that produces only one file in the end, using a list of functions (with arguments)

    # Remove all intermittent results just in case!!!
    # extract surfaces
    # apply_single_transformation_to_all(path, input_base='h_case', version='', start=start_, end=end_, ext='_',
    #                                    ext_type='PolyData', function_='extract_surface')
    #
    # # scale the meshes
    # apply_single_transformation_to_all(path, input_base='h_case', version='_', start=start_, end=end_, ext='',
    #                                    ext_type='UG', function_='scale', args='((.001, .001, .001))')
    # center the meshes
    # apply_single_transformation_to_all(path, input_base='h_case', version='', start=start_, end=end_, ext='',
    #                                    ext_type='UG', function_='translate_to_center')
    # clean (if possible) poorly built meshes
    # apply_single_transformation_to_all(path, input_base='h_case', version='_', start=start_, end=end_, ext='',
    #                                    ext_type='UG', function_='clean_polydata', args='(1e-6, True)')
    # align the meshes
    apply_function_to_all(path, input_base='h_case', version='', start=start_, end=end_, ext='algn',
                          ext_type='UG', function_='alignment', args='labels=(7, 8), ')
    # split chambers
    apply_function_to_all(path, 'h_case', version='algn',  start=start_, end=end_, ext='_surface_full',
                          function_='split_chambers', args='case, return_elements=True')
    # -----Slice extraction-----
    # create slices
    # apply_function_to_all(path, input_base='h_case', version='', start=start_, end=end_, ext='plax',
    #                       function_='create_plax_slices', args='')
    # # save slices as png
    # apply_single_transformation_to_all(path, input_base='h_case', version='plax', start=start_, end=end_, ext=None,
    #                                    function_='write_png', args='()')
    # decimate_heart
    # apply_function_to_all(path, 'case_', '_pd_centered', start=start_, end=end_, ext='_delete',
    # function_='decimate_heart',
    # args='case')
    pass


if __name__ == '__main__':

    # deformetrica_data_path = os.path.join('/home', 'mat', 'Deformetrica')
    # def_relevant_files = glob.glob(os.path.join(deformetrica_data_path,
    #                                             'deterministic_atlas_ct',
    #                                             'output_tmp10_def10_surf',
    #                                             'DeterministicAtlas__flow__heart__subject_sub??__tp_?.vtk'))
    # def_relevant_files.sort()

    absolute_data_path = os.path.join('/home', 'mat', 'Python', 'data', 'h_case')
    relevant_files = glob.glob(os.path.join(absolute_data_path, 'h_case*.vtk'))
    relevant_files.sort()

    tetra_data_path = os.path.join('/home', 'mat', 'Deformetrica', 'deterministic_atlas_ct', 'gmsh', 'tetra_template')
    tetra_files = glob.glob(os.path.join(tetra_data_path, '*tetra*'))
    tetra_files.sort()

# -----Main Pipeline
#     h_case_pipeline(path=absolute_data_path, start_=1, end_=19)
# -----------------

# ----Building models for electrophysiological simulations
#     models = []
#     for element in tetra_files:
#         print(element)
#         print(os.path.basename(element))
#         element_name = os.path.basename(element).split('_')[0]
#         print(element_name)
#         model = Model(element)
#         element_tag = [i+1 for i, elem in enumerate(model.list_of_elements) if elem == element_name][0]
#         print('Element name: {}, element tag: {}'.format(element_name, element_tag))
#         model.build_tag(label=element_tag)
#         model.change_tag_label()
#         models.append(model)
#
#     final_model = models.pop(0)
#     for model_to_merge in models:
#         final_model.mesh = merge_elements(final_model.mesh, model_to_merge.mesh)
#     final_model.tetrahedralize()
#     final_model.write_vtk(postscript='merged', type_='UG')
# -----------------------------------------------------------------------

# -----Testing single function
#     relevant_files = [x for x in relevant_files if 'plax' in x]
#     print(relevant_files)
#     model = Model(filename=relevant_files[0])
#     for lv in relevant_files:
#         model = Model(lv)
#         get_lowest_septal_point(model)
# ----------------------------

# -----Frangi dataset------------------------------------------------------------------------------------------
# labels_and_ranges = ({'label': 1, 'range': (0, 472726)}, {'label': 2, 'range': (472726, 688801)})
# model.mesh = assign_tags(model.mesh, labels_and_ranges)
# model.write_vtk(postscript='tagged', type_='UG')

# model = Model('h_case/h_case01_surface_la.vtk')  # Relative, path to the file
# iterations = 200
# model.smooth_laplacian(iterations)
# # model.normals()
# model.write_vtk(postscript='_smooth_' + str(iterations))
# decimate_chambers()

# model.smooth_window(100, 0.0001)
# model.normals()
# model.write_vtk(postscript='smooth')

# model = Model("/home/mat/Deformetrica/deterministic_atlas_ct/Temp/DeterministicAtlas__EstimatedParameters__Template_LV.vtk")

# # This is how we define modes. They come in pairs
# # 'name of mode': scale
# # names are from mode_01 to mode_41
# # scale can be positive or negative, 0 means no deformation in that direction
# # modes can be added or removed
# # modes = {'mode_01': -2.5, 'mode_02': -2.4, 'mode_03': 2.2, 'mode_07': 1.2}
# # model.apply_modes(modes)
# model.visualize_mesh()
# model.write_vtk(postscript='warped')
# ------------------------------------------------------------------------------------------------------------
