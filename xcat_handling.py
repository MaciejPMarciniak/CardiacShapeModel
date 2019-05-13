import numpy as np
import matplotlib.pyplot as plt


class XCat:

    image_shape = [1800, 400, 600]

    def __init__(self, filename=''):
        self.filename = filename
        self.image = None
        self.slice = None

    def _read_file(self):

        image = np.fromfile(self.filename, dtype='>i2', sep='')
        self.image = image.reshape(self.image_shape)  # The order of voxel numbers is reversed

    def _get_slice(self, slice_number=0, view=''):

        if self.image is None:
            self._read_file()

        if view == 'front':
            slice_ = self.image[:, slice_number, :]
        elif view == 'side':
            slice_ = self.image[:, :, slice_number]
        else:
            slice_ = self.image[slice_number, :, :]
        return slice_

    def save_single_slice(self, slice_number=300, view=''):

        csv_file = self.filename.split('.')[0] + '.csv'
        np.savetxt(csv_file, self._get_slice(slice_number, view), delimiter=',')

    def show_single_slice(self, slice_number=200, view=''):

        plt.imshow(self._get_slice(slice_number, view), cmap='Greys_r')
        plt.show()


if __name__ == '__main__':

    young = '/home/mat/Desktop/models/10yr_male_600_400_1800_16bit_BE.raw'
    adult = '/home/mat/Desktop/models/adult_male_600_400_1800_16bit_BE.raw'

    xcat = XCat(young)
    xcat.save_single_slice(700)
    xcat.show_single_slice(900)
