import os
import glob
from shutil import copyfile, move
import subprocess
from VTK_background import Heart, merge_elements
from pathlib import Path


class ShootingSupport:

    def __init__(self, main_path, models_path, geo_path, temp_path, output_path, k_model, template=False):
        self.main_path = main_path
        self.models_path = models_path
        self.temp_path = temp_path
        self.geo_path = geo_path
        self.output_path = output_path
        self.k_model = k_model
        self.template = template

    def clean(self):
        [f.unlink() for f in Path(self.temp_path).glob("*") if f.is_file()]

    def copy_model(self):
        self.clean()
        print(self.models_path)
        if not self.template:
            model_files = glob.glob(os.path.join(self.models_path, 'Shooting_{}_*'.format(self.k_model)))
        else:
            model_files = glob.glob(os.path.join(models_path, '*Template*'))
        print(model_files)
        model_files = [mf for mf in model_files if 'ControlPoints' not in mf and'Momenta' not in mf]
        for mf in model_files:
            copyfile(mf, os.path.join(self.temp_path, os.path.basename(mf)))

    def modify_geo_files(self):
        geo_files = glob.glob(os.path.join(self.geo_path, '*.geo'))
        geo_files.sort()
        model_files = glob.glob(os.path.join(self.temp_path, '*'))
        model_files.sort()

        for gf, mf in zip(geo_files, model_files):
            with open(gf, 'r') as file:
                data = file.readlines()
            data[1] = ('Merge "' + mf + '";\n')
            with open(gf, 'w') as file:
                file.writelines(data)

    def run_tetrahedralization(self):
        subprocess.call(os.path.join(self.main_path, './meshing.sh'))

    def tag_and_merge(self):
        models = []
        tetra_files = [os.path.join(self.main_path, 'tetra', el + '_tetra.vtk') for el in Heart.list_of_elements]
        for element in tetra_files:
            print(os.path.basename(element))
            element_name = os.path.basename(element).split('_')[0]
            print(element_name)
            model = Heart(element)
            element_tag = [i + 1 for i, elem in enumerate(model.list_of_elements) if elem == element_name][0]
            print('Element name: {}, element tag: {}'.format(element_name, element_tag))
            model.build_tag(label=element_tag)
            model.change_tag_label()
            models.append(model)

        final_model = models.pop(0)
        for model_to_merge in models:
            final_model.mesh = merge_elements(final_model.mesh, model_to_merge.mesh)
        final_model.tetrahedralize()
        final_model.write_vtk(postscript='merged', type_='UG')
        if not self.template:
            os.rename(os.path.join(self.main_path, 'tetra', 'LV_tetramerged.vtk'),
                      os.path.join(self.main_path, 'tetra', 'Full_Heart_{}.vtk'.format(self.k_model)))
            move(os.path.join(self.main_path, 'tetra', 'Full_Heart_{}.vtk'.format(self.k_model)),
                 os.path.join(self.output_path, 'Full_Heart_{}.vtk'.format(self.k_model)))
        else:
            os.rename(os.path.join(self.main_path, 'tetra', 'LV_tetramerged.vtk'),
                      os.path.join(self.main_path, 'tetra', 'Full_Template.vtk'))
            move(os.path.join(self.main_path, 'tetra', 'Full_Template.vtk'),
                 os.path.join(self.output_path, 'Full_Template.vtk'))

    def pipeline_surf_2_tetra_mesh(self):
        self.copy_model()
        self.modify_geo_files()
        self.run_tetrahedralization()
        self.tag_and_merge()


if __name__ == '__main__':

    # TODO: Make sure that the element files of vtk and geo are sorted in the same way in 'modify_geo_files'

    for i in range(0, 1):
        output_path = '/media/mat/D6B7-122E/Final_models_{}'.format(str(i).zfill(2))
        models_path = ('/home/mat/Deformetrica/deterministic_atlas_ct/'
                       'output_shooting/extreme_shapes').format(i)
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        for j in range(0, 18):

            sup = ShootingSupport(main_path='/home/mat/Deformetrica/deterministic_atlas_ct/gmsh',
                                  models_path=models_path,
                                  geo_path='/home/mat/Deformetrica/deterministic_atlas_ct/gmsh/geofiles',
                                  temp_path='/home/mat/Deformetrica/deterministic_atlas_ct/gmsh/current_model',
                                  output_path=output_path,
                                  k_model=j,
                                  template=False)
            sup.pipeline_surf_2_tetra_mesh()
