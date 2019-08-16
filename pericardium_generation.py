from VTK_background import Model, merge_elements


def create_pericarium(path_to_vtk_file=''):

    heart = Model(path_to_vtk_file)
    heart.extract_surface()
    heart.get_external_surface()
    heart.unstructured_grid_to_poly_data()

    peri = Model(path_to_vtk_file)
    peri.extract_surface()
    peri.implicit_modeller(1)
    peri.contouring()
    peri.get_external_surface()
    peri.smooth_window(30, 0.05)
    peri.unstructured_grid_to_poly_data()
    peri.mesh = merge_elements(heart.mesh, peri.mesh)
    peri.write_vtk()


if __name__ == '__main__':
    vtk_file = "/home/mat/Deformetrica/deterministic_atlas_ct/Main_modes_export/Template_tetra.vtk"
    create_pericarium(vtk_file)
