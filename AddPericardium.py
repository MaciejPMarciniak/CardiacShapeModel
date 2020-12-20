from Mesh import Model


def create_pericarium(path_to_vtk_file=''):

    heart = Model(path_to_vtk_file)
    heart.extract_surface()
    heart.get_external_surface()
    heart.unstructured_grid_to_poly_data()
    heart.write_vtk(postscript='_mesh_surf')

    peri = Model(path_to_vtk_file)
    peri.extract_surface()
    peri.implicit_modeller(1)
    peri.contouring()
    peri.get_external_surface()
    peri.unstructured_grid_to_poly_data()
    peri.write_vtk(postscript='_peri_surf')


if __name__ == '__main__':
    vtk_file = "C:\Data\DataGeneration\Extreme3\mode1neg.vtk"
    create_pericarium(vtk_file)
