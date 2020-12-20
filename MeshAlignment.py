import numpy as np
from vtk.util.numpy_support import vtk_to_numpy


def calculate_rotation(reference_vector, target_vector):
    """
    :param reference_vector: Vector with referential direction. The rotation matrix will align the target_vector's
    direction to this one.
    :param target_vector: Target mesh direction vector.
    :return: 3x3 rotation matrix (rot), where [rot @ target_vector = reference_vector]
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
    """
    :param _model: Mesh.Model object
    :param _labels: labels of elements to calculate the centers of
    :return: Centers of mass of the mesh elements
    """
    centers = []
    for lab in _labels:
        centers.append(_model.get_center(_model.threshold(lab, lab)))
    print('centers: {}'.format(centers))
    return centers


def get_translation_vector(target_markers, reference_markers):
    """
    :param target_markers: Two 3D points of the target mesh corresponding to the reference markers
    :param reference_markers: Two 3D points used for positional alignment
    :return: translation vector to bring the target mesh to the reference mesh position
    """
    target_center = np.mean([target_markers[0], target_markers[1]], axis=0)
    reference_center = np.mean([reference_markers[0], reference_markers[1]], axis=0)
    return reference_center - target_center


def get_vector_alignment_rotation_matrix(target_markers, reference_markers):
    """
        :param target_markers: Two 3D points of the target mesh corresponding to the reference markers
        :param reference_markers: Two 3D points used for positional alignment
        :return: rotation matrix to bring the target mesh to the reference mesh orientation
    """
    target_center = np.mean([target_markers[0], target_markers[1]], axis=0)
    reference_center = np.mean([reference_markers[0], reference_markers[1]], axis=0)
    return calculate_rotation(reference_markers[2] - reference_center, target_markers[2] - target_center)


def get_plane_alignment_rotation_matrix(target_markers, reference_markers):
    """
    :param target_markers: Three 3D points of the target mesh corresponding to the reference markers
    :param reference_markers: Three 3D points used for positional alignment
    :return: rotation matrix to bring the target mesh to the reference mesh orientation
    """
    target_center = np.mean([target_markers[0], target_markers[1]], axis=0)
    reference_center = np.mean([reference_markers[0], reference_markers[1]], axis=0)
    assert np.all(target_center - reference_center < 0.001), 'The models are not position-aligned'

    target_plane_normal = calculate_plane_normal(target_markers[2], target_center, target_markers[1])
    reference_plane_normal = calculate_plane_normal(reference_markers[2], reference_center, reference_markers[1])
    return calculate_rotation(reference_plane_normal, target_plane_normal)


def get_lowest_septal_point(_model):
    """
    :param _model: Mesh.Model object
    :return: Find the reference point for rotation and translation (the point with maximum distance to LV and RV valves)
    """
    lv = _model.threshold(1, 1)
    rv = _model.threshold(2, 2)
    lv_points = set(tuple(map(tuple, vtk_to_numpy(lv.GetOutput().GetPoints().GetData()))))
    rv_points = set(tuple(map(tuple, vtk_to_numpy(rv.GetOutput().GetPoints().GetData()))))
    common_points = lv_points.intersection(rv_points)

    valve_centers = get_centers(_model, (7, 8))
    center = np.mean([valve_centers[0], valve_centers[1]], axis=0)
    norms = np.array([[x, np.linalg.norm((center - x))] for x in common_points])
    lowest_septal_point = np.array(norms[norms[:, 1] == np.max(norms[:, 1])][0, 0])
    return lowest_septal_point


def alignment(target_model, reference_model, labels=(7, 8)):
    """
    Rigidly align the two meshes by bringing the target model to the same space as the reference model. The landmarks
    are the centers of mass of the mitral and tricuspid valves and the point that is most distant from them within the
    interventricular septum.
    :param target_model: model to be aligned
    :param reference_model: model to which all the models are aligned within the cohort
    :param labels: tags of the mitral and tricuspid valves in the meshes
    :return: aligned target model
    """
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
