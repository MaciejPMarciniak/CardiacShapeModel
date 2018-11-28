import vtk
import numpy as np
import os
import glob
import fnmatch


class Heart:

    def __init__(self, filename='h_case06.vtk', to_polydata=False):

        self.filename, self.input_type = filename.split('.')
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

        self.scalar_range = [1.0, 17.0]  # Added for particular case of CT meshes
        print('Corrected scalar range: {}'.format(self.scalar_range))
        self.center_of_heart = self.set_center(self.mesh)

    @staticmethod
    def set_center(_mesh):
        centerofmass = vtk.vtkCenterOfMass()
        centerofmass.SetInputData(_mesh.GetOutput())
        centerofmass.Update()
        return np.array(centerofmass.GetCenter())

    @staticmethod
    def unstructured_grid_to_poly_data(_algorithm):
        geometry_filter = vtk.vtkExtractGeometry()
        geometry_filter.SetInputConnection(_algorithm.GetOutputPort())
        geometry_filter.Update()
        return geometry_filter

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
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # Create the Renderer
        renderer = vtk.vtkRenderer()
        renderer.ResetCameraClippingRange()
        renderer.AddActor(actor)  # More actors can be added
        renderer.SetBackground(1, 1, 1)  # Set background to white

        # Create the RendererWindow
        renderer_window = vtk.vtkRenderWindow()
        renderer_window.AddRenderer(renderer)

        # Display the mesh
        # noinspection PyArgumentList
        if display:
            interactor = vtk.vtkRenderWindowInteractor(actually='', one_of_its_subclasses='')
            interactor.SetRenderWindow(renderer_window)
            interactor.Initialize()
            interactor.Start()
        else:
            return renderer_window

    # -----3D rigid transformations---------------------------------------------------------------------------
    def translate_to_center(self):
        # vtkTransform.SetMatrix - enables for applying 4x4 transformation matrix to the meshes
        translate = vtk.vtkTransform()
        translate.Translate(-self.center_of_heart[0], -self.center_of_heart[1], -self.center_of_heart[2])
        translate.Update()
        transformer = vtk.vtkTransformFilter()
        transformer.SetInputConnection(self.mesh.GetOutputPort())
        transformer.SetTransform(translate)
        transformer.Update()
        self.mesh = transformer
        self.center_of_heart = self.set_center(self.mesh)
        print('Center:', self.center_of_heart)

    def rotate(self, alpha=0, beta=0, gamma=0):
        rotate = vtk.vtkTransform()
        rotate.Identity()
        rotate.RotateZ(-45)
        transformer = vtk.vtkTransformFilter()
        transformer.SetInputConnection(self.mesh.GetOutputPort())
        transformer.SetTransform(rotate)
        transformer.Update()
        self.mesh = transformer
        self.center_of_heart = self.set_center(self.mesh)

    def scale(self, factor=(0.001, 0.001, 0.001)):
        scale = vtk.vtkTransform()
        scale.Scale(factor[0], factor[1], factor[2])
        transformer = vtk.vtkTransformFilter()
        transformer.SetInputConnection(self.mesh.GetOutputPort())
        transformer.SetTransform(scale)
        transformer.Update()
        self.mesh = transformer
        self.center_of_heart = self.set_center(self.mesh)

    # -----Mesh manipulation----------------------------------------------------------------------------------
    def downsample(self, tolerance=0.05):
        downsample = vtk.vtkCleanPolyData()
        downsample.SetInputConnection(self.mesh.GetOutputPort())
        downsample.SetTolerance(tolerance)
        downsample.Update()
        self.mesh = downsample

    def threshold(self, low=0, high=100):
        threshold = vtk.vtkThreshold()
        threshold.SetInputConnection(self.mesh.GetOutputPort())
        threshold.ThresholdBetween(low, high)
        threshold.Update()
        return threshold

    def smooth_laplacian(self, number_of_iterations=50):
        smooth = vtk.vtkSmoothPolyDataFilter()
        smooth.SetInputConnection(self.mesh.GetOutputPort())
        smooth.SetNumberOfIterations(number_of_iterations)
        smooth.FeatureEdgeSmoothingOff()
        smooth.BoundarySmoothingOn()
        smooth.Update()
        self.mesh = smooth

    def smooth_window(self, number_of_iterations=15, pass_band=0.5):
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

    def delaunay3d(self):
        delaunay3d = vtk.vtkDelaunay3D()
        delaunay3d.SetInputConnection(self.mesh.GetOutputPort())
        delaunay3d.Update()
        self.mesh = delaunay3d

    def delaunay2d(self):
        delaunay2d = vtk.vtkDelaunay2D()
        delaunay2d.SetInputConnection(self.mesh.GetOutputPort())
        delaunay2d.Update()
        self.mesh = delaunay2d

    def normals(self):
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(self.mesh.GetOutputPort())
        normals.FlipNormalsOn()
        normals.Update()
        self.mesh = normals

    def decimation(self, reduction=50):
        decimation = vtk.vtkQuadricDecimation()
        decimation.SetInputConnection(self.mesh.GetOutputPort())
        decimation.VolumePreservationOn()
        decimation.SetTargetReduction(reduction / 100)  # percent of removed triangles
        decimation.Update()
        self.mesh = decimation

    def subdivision(self, number_of_subdivisions=3):
        self.normals()
        subdivision = vtk.vtkLinearSubdivisionFilter()
        subdivision.SetNumberOfSubdivisions(number_of_subdivisions)
        subdivision.SetInputConnection(self.mesh.GetOutputPort())
        subdivision.Update()
        self.mesh = subdivision
        self.visualize_mesh(True)

    def fill_holes(self, hole_size=10.0):
        filling_filter = vtk.vtkFillHolesFilter()
        filling_filter.SetInputConnection(self.mesh.GetOutputPort())
        filling_filter.SetHoleSize(hole_size)
        filling_filter.Update()
        self.mesh = filling_filter

    def apply_modes(self, modes_with_scales):
        for mode, scale in modes_with_scales.items():
            print('Applying ' + mode + ' multiplied by ' + str(scale))
            self.mesh.GetOutput().GetPointData().SetActiveVectors(mode)
            warp_vector = vtk.vtkWarpVector()
            warp_vector.SetInputConnection(self.mesh.GetOutputPort())
            warp_vector.SetScaleFactor(scale)
            warp_vector.Update()
            self.mesh = warp_vector

    def extract_surface(self):
        # Get surface of the mesh
        surface_filter = vtk.vtkDataSetSurfaceFilter()
        surface_filter.SetInputData(self.mesh.GetOutput())
        surface_filter.Update()
        self.mesh = surface_filter

    @staticmethod
    def get_external_surface(_mesh, external=True):
        _center = np.zeros(3)
        _bounds = np.zeros(6)
        _ray_start = np.zeros(3)
        cell_id = vtk.mutable(-1)
        xyz = np.zeros(3)
        pcoords = np.zeros(3)
        t = vtk.mutable(0)
        sub_id = vtk.mutable(0)
        if external:
            surf = 1.1
        else:
            surf = -1.1

        _mesh.GetOutput().GetCenter(_center)
        _mesh.GetOutput().GetPoints().GetBounds(_bounds)
        for j in range(3):
            _ray_start[j] = _bounds[2 * j + 1] * surf

        cell_locator = vtk.vtkCellLocator()
        cell_locator.SetDataSet(_mesh.GetOutput())
        cell_locator.BuildLocator()
        cell_locator.IntersectWithLine(_ray_start, _center, 0.0001, t, xyz, pcoords, sub_id, cell_id)
        print('ID of the cell on the outer surface: {}'.format(cell_id))

        connectivity_filter = vtk.vtkConnectivityFilter()
        connectivity_filter.SetInputConnection(_mesh.GetOutputPort())
        connectivity_filter.SetExtractionModeToCellSeededRegions()
        connectivity_filter.InitializeSeedList()
        connectivity_filter.AddSeed(cell_id)
        connectivity_filter.Update()
        return connectivity_filter

    # -----MeshInformation------------------------------------------------------------------------------------
    def print_numbers(self):
        _mesh = self.mesh.GetOutput()
        print('Number of verices: {}'.format(_mesh.GetNumberOfVerts()))
        print('Number of lines: {}'.format(_mesh.GetNumberOfLines()))
        print('Number of strips: {}'.format(_mesh.GetNumberOfStrips()))
        print('Number of polys: {}'.format(_mesh.GetNumberOfPolys()))
        print('Number of cells: {}'.format(_mesh.GetNumberOfCells()))
        print('Number of points: {}'.format(_mesh.GetNumberOfPoints()))

    def get_volume(self):
        mass = vtk.vtkMassProperties()
        mass.SetInputConnection(self.mesh.GetOutputPort())
        return mass.GetVolume()

    # -----InputOutput----------------------------------------------------------------------------------------
    def read_vtk(self, to_polydata=False):
        # Read the source file.
        reader = vtk.vtkDataReader()
        reader.SetFileName(self.filename + '.' + self.input_type)
        reader.Update()

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
        reader.SetFileName(self.filename + '.' + self.input_type)
        reader.Update()
        scalar_range = reader.GetOutput().GetScalarRange()
        return reader, scalar_range

    def read_obj(self):
        reader = vtk.vtkOBJReader()
        reader.SetFileName(self.filename + '.' + self.input_type)
        reader.Update()
        scalar_range = reader.GetOutput().GetScalarRange()
        return reader, scalar_range

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

    def write_obj(self):
        output_filename = self.filename
        render_window = self.visualize_mesh(False)

        print('Saving PolyData in the OBJ file...')
        obj_writer = vtk.vtkOBJExporter()
        obj_writer.SetRenderWindow(render_window)
        obj_writer.SetFilePrefix(output_filename)
        obj_writer.Write()
        print('{} written succesfully'.format(output_filename + '.obj'))

    def write_vtk(self, postscript='_new', type_='PolyData'):
        output_filename = self.filename + postscript + '.vtk'
        writer = None
        if type_ == 'PolyData':
            print('Saving PolyData...')
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


# TODO: Make all of the functions and parameters below into a class!!!
# -----MightBeUseful------------------------------------------------------------------------------------------
def get_valve_surfaces(_model):

    valve = _model.threshold(14, 14)  # mitral valve
    valve_ext_left = _model.unstructured_grid_to_poly_data(_model.get_external_surface(valve))
    valve_int_left = _model.unstructured_grid_to_poly_data(_model.get_external_surface(valve, False))
    valve_int_left = change_elem_tag(valve_int_left, 25)
    valve = _model.threshold(15, 15)  # triscupid valve
    valve_ext_right = _model.unstructured_grid_to_poly_data(_model.get_external_surface(valve))
    valve_int_right = _model.unstructured_grid_to_poly_data(_model.get_external_surface(valve, False))
    valve_int_right = change_elem_tag(valve_int_right, 26)
    lower_part = _model.unstructured_grid_to_poly_data(_model.threshold(1, 13))
    upper_part = _model.unstructured_grid_to_poly_data(_model.threshold(16, 24))

    app = vtk.vtkAppendPolyData()
    app.AddInputConnection(lower_part.GetOutputPort())
    app.AddInputConnection(upper_part.GetOutputPort())
    app.AddInputConnection(valve_ext_left.GetOutputPort())
    app.AddInputConnection(valve_int_left.GetOutputPort())
    app.AddInputConnection(valve_ext_right.GetOutputPort())
    app.AddInputConnection(valve_int_right.GetOutputPort())
    app.Update()
    _model.mesh = app
    return _model


# -----ApplyToCohort------------------------------------------------------------------------------------------
def apply_single_transformation_to_all(input_base, version, start=1, end=20, ext='_new', function_=None, args='()'):
    if function_ is not None:
        for case_no in range(start, end+1):
            case = input_base + '/' + input_base + str(case_no).zfill(2) + version + '.vtk'
            print(case)
            single_model = Heart(case)
            exec('model.' + function_ + args)
            single_model.write_vtk(ext)


def apply_function_to_all(input_base, version, start=1, end=20, ext='_new', function_=None):
    if function_ is not None:
        for case_no in range(start, end+1):
            case = input_base + '/' + input_base + str(case_no).zfill(2) + version + '.vtk'
            print(case)
            single_model = Heart(case)
            single_model = function_(single_model, case)
            single_model.write_vtk(ext)
# ------------------------------------------------------------------------------------------------------------


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
    FINAL_POINTS = 10000
    built_points = 0
    reduction = (1 - FINAL_POINTS / all_points) * 100
    chambers = ['lv', 'rv', 'ra', 'la']
    case_split = case.split('.')
    sum_chamber_points = 0
    for chamber in chambers:
        print('decimating ' + chamber)
        _model = Heart(case_split[0].strip('ful') + chamber + '.' + case_split[1])
        _model.smooth_laplacian(500)
        chamber_points = _model.mesh.GetOutput().GetNumberOfPoints()
        expected_number_of_points = np.round((1 - reduction / 100) * chamber_points)
        built_points += expected_number_of_points
        if chamber == 'la':
            expected_number_of_points += FINAL_POINTS - built_points
        reduction = hill_climb_with_magnitude(initial_x=reduction, eval_function=eval_decimation,
                                              goal=expected_number_of_points, step=0.001, this_model=_model)
        _model.decimation(reduction=reduction)
        sum_chamber_points += _model.mesh.GetOutput().GetNumberOfPoints()
        _model.normals()
        _model.write_vtk(postscript='_decimated')
        # exit('checking decimation in decimate_chambers!')
    if sum_chamber_points != FINAL_POINTS:
        exit('Chamber_points: {}, but should be {}'.format(sum_chamber_points, FINAL_POINTS))

    return full_model
# ------------------------------------------------------------------------------------------------------------


def change_elem_tag(_mesh, label):
    size = _mesh.GetOutput().GetAttributes(1).GetArray(0).GetSize()
    for i in range(size):
        _mesh.GetOutput().GetAttributes(1).GetArray(0).SetTuple(i, (float(label),))
    return _mesh


def split_chambers(_model, return_as_surface=True):
    _model.translate_to_center()

    surfaces = []
    for i in range(1, int(_model.scalar_range[1]) + 1):

        x = _model.threshold(i, i)
        surfaces.append(x)

    full_model_appended = vtk.vtkAppendFilter()
    for surf in surfaces:
        full_model_appended.AddInputConnection(surf.GetOutputPort())
        # _
    full_model_appended.Update()
    _model.mesh = full_model_appended
    if return_as_surface:
        _model.extract_surface()
        _model.write_vtk(postscript='surf')
        _model.write_obj()
    else:
        _model.write_vtk(postscript='tetra')
    return _model


def split_and_combine_chambers(_model, case=None):
    # Thresholds related to arteries commented, as well as saving particular meshes
    _model.translate_to_center()
    _mesh = _model.mesh
    # labels of chambers
    lv = [1, 7, 9]
    rv = [2, 8, 10]
    la = [3, 11, 12, 13, 14, 15]
    ra = [4, 16, 17]
    aorta = [5]
    pa = [6]
    
    surfaces = {'surface_lv': [], 'surface_rv': [], 'surface_la': [], 'surface_ra': [],
                'surface_aorta': [], 'surface_pa': []}

    for i in range(1, int(_model.scalar_range[1])+1):
        surface = _model.threshold(i, i)
        if i in lv:
            surfaces['surface_lv'].append(surface)
        elif i in rv:
            surfaces['surface_rv'].append(surface)
        elif i in la:
            surfaces['surface_la'].append(surface)
        elif i in ra:
            surfaces['surface_ra'].append(surface)
        elif i in aorta:
            surfaces['surface_aorta'].append(surface)
        elif i in pa:
            surfaces['surface_pa'].append(surface)
        else:
            print('Label {} doesn\'t belong to any chamber of interest'.format(i))

    full_model_appended = vtk.vtkAppendPolyData()
    for surf, part in surfaces.items():
        app = vtk.vtkAppendFilter()
        for element in part:
            app.AddInputConnection(element.GetOutputPort())
        app.Update()
        _model.mesh = app
        _model.extract_surface()
        full_model_appended.AddInputConnection(_model.mesh.GetOutputPort())
        # _model.visualize_mesh()
        # _model.write_vtk(postscript=surf)
    full_model_appended.Update()
    _model.mesh = full_model_appended

    return _model


def change_downloaded_files_names(path='h_case_', key='surfmesh', ext='vtk'):
    files = glob.glob(os.path.join(path, '*'+key+'*'+ext))
    for i, old_file in enumerate(files):
        new_file = old_file.split('.')[0]
        print(new_file)
        os.rename(old_file, new_file+'.'+ext)
    print(files)


def h_case_pipeline(start_=1, end_=20):
    # TODO: Problem with this solution is that the files are read and written at every step. Should become
    # TODO: a pipeline, that produces only one file in the end, using a list of functions (with arguments)

    # extract surfaces
    # apply_single_transformation_to_all('case_', version='', start=start_, end=end_, ext='_pd',
    #                                    function_='extract_surface')
    # # center the meshes
    apply_single_transformation_to_all('case_', version='_pd', start=start_, end=end_, ext='_centered',
                                       function_='translate_to_center')
    # scale the meshes
    apply_single_transformation_to_all('case_', version='_pd_centered', start=start_, end=end_, ext='',
                                       function_='scale', args='()')
    # extract valves
    # apply_function_to_all('h_case', '',  start=start_, end=end_, ext='_sepvalves', function_=get_valve_surfaces)
    # split chambers
    apply_function_to_all('h_case_', '_surface_pd_centered',  start=start_, end=end_,ext='_surface_full', function_=split_chambers)
    # decimate_heart
    apply_function_to_all('case_', '_pd_centered', start=start_, end=end_, ext='_delete', function_=decimate_heart)
    # write as mha
    # apply_single_transformation_to_all('h_case', '_surface_lv_decimated', start=start_, end=end_, function_='write_mha()')
    pass


if __name__ == '__main__':
    absolute_data_path = os.path.join('/home', 'mat', 'Python', 'data', )
    #    h_case_pipeline(start_=18, end_=18)
    relevant_files = glob.glob(os.path.join(absolute_data_path, 'case_', 'tetra_meshes', 'h_case??surf_doubled.vtk'))
    relevant_files.sort()
    print(relevant_files)
    for i, shape in enumerate(relevant_files):

        model = Heart(shape)
        # model.write_vtk(postscript='surf_doubled')
        model.write_stl()
        model.write_obj()


    # model = Heart('/home/mat/Deformetrica/deterministic_atlas_ct/output_tmp_10_def_10/DeterministicAtlas__flow__heart__subject_sub01__tp_10.vtk')
    #     model.write_stl()

    # model = Heart('h_case_/h_case_01_surface_pd_centered.vtk')  # Relative, path to the file

# # --Attempt to translate unstructured grid into mha image format----------
#     inval = 255
#     outval = 0
#     spacing = np.ones(3) * 0.05
#     dim = np.zeros(3, dtype='i8')
#     origin = np.zeros(3)
#     bounds = np.zeros(6)
#     pd = model.mesh.GetOutput().GetBounds(bounds)
#
#     for i in range(3):
#         dim[i] = np.ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i])
#     for i in range(3):
#         origin[i] = bounds[i*2] + spacing[i] / 2
#
#     print(dim)
#     print(origin)
#
#     whiteImage = vtk.vtkImageData()
#     whiteImage.SetSpacing(spacing)
#     whiteImage.SetDimensions(dim)
#     whiteImage.SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1)
#     whiteImage.SetOrigin(origin)
#     whiteImage.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)
#     pti = vtk.vtkPolyDataToImageStencil()
#     pti.SetInputData(model.mesh.GetOutput())
#     pti.Update()
#     stencil = vtk.vtkImageStencil()
#     stencil.SetInputData(whiteImage)
#     stencil.SetStencilData(pti.GetOutput())
#     stencil.ReverseStencilOff()
#     stencil.SetBackgroundValue(1)
#     stencil.Update()
#
#     writer = vtk.vtkMetaImageWriter()
#     writer.SetFileName("SphereVolume.mhd");
#     writer.SetInputData(stencil.GetOutput());
#     writer.Write()
#  -----version 2---------------

#     count = whiteImage.GetNumberOfPoints()
#     print(count)
#     for i in range(count):
#         print('Points added: {}%'.format(np.round(i/count, 2))) if i % 100000 == 0 else None
#         whiteImage.GetPointData().GetScalars().SetTuple(i, (float(inval),))
#
#     print('PolyData to ImageStencil...')
#     pol2stenc = vtk.vtkPolyDataToImageStencil()
#     pol2stenc.SetInputData(pd)
#     pol2stenc.SetOutputOrigin(origin)
#     pol2stenc.SetOutputSpacing(spacing)
#     pol2stenc.SetOutputWholeExtent(whiteImage.GetExtent())
#     pol2stenc.Update()
#
#     print('Image Stencil...')
#     imgstenc = vtk.vtkImageStencil()
#     imgstenc.SetInputData(whiteImage)
#     imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
#     imgstenc.ReverseStencilOff()
#     imgstenc.SetBackgroundValue(outval)
#     imgstenc.Update()
#
#     print('Writing the file...')
#     writer = vtk.vtkMetaImageWriter()
#     writer.SetFileName("SphereVolume.mha");
#     writer.SetInputData(imgstenc.GetOutput());
#     writer.Write()
# ------------------------------------------------------------------------------------------------------------
# model = Heart('h_case/h_case01_surface_la.vtk')  # Relative, path to the file
# iterations = 200
# model.smooth_laplacian(iterations)
# # model.normals()
# model.write_vtk(postscript='_smooth_' + str(iterations))
# decimate_chambers()

# model.smooth_window(100, 0.0001)
# model.normals()
# model.write_vtk(postscript='smooth')

# model.visualize_mesh()  # Mean shape. Close the mesh window to move on
# model.get_external_surface()
# model.write_vtk()
# model.write_obj()
# This is how we define modes. They come in pairs
# 'name of mode': scale
# names are from mode_01 to mode_41
# scale can be positive or negative, 0 means no deformation in that direction
# modes can be added or removed
# modes = {'mode_01': -2.5, 'mode_02': -2.4, 'mode_03': 2.2, 'mode_07': 1.2}
# model.apply_modes(modes)
# model.visualize_mesh()

# ------------------------------------------------------------------------------------------------------------
