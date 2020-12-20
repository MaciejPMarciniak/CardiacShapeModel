import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from MeshAlignment import get_centers, calculate_plane_normal


def get_apical_landmarks(_model, valve_labels):
    lv = _model.threshold(1, 1)
    lv_points = vtk_to_numpy(lv.GetOutput().GetPoints().GetData())

    valve_centers = get_centers(_model, valve_labels)
    apex_id = np.argmax([np.linalg.norm(valve_centers[0] - x) for x in lv_points])
    apex = lv_points[apex_id]

    normal = calculate_plane_normal(*valve_centers, apex)
    origin = np.mean(np.array((apex, *valve_centers)), axis=0)

    return origin, normal, (*valve_centers, apex)


# ---Parasternal long axis view------------------------------------------------------------------------------
def create_plax_slices(_model):
    origin, normal, landmarks = get_apical_landmarks(_model, (7, 9))
    _model.slice_extraction(origin, normal)
    _model.align_slice(landmarks[2], landmarks[1], landmarks[0])
    _model.rotate(gamma=-90)
    return _model
# ---END-Parasternal long axis view--------------------------------------------------------------------------


# ---Four chamber view---------------------------------------------------------------------------------------
def create_4ch_slices(_model):
    origin, normal, landmarks = get_apical_landmarks(_model, (7, 8))
    _model.slice_extraction(origin, normal)
    _model.align_slice(landmarks[2], landmarks[1], landmarks[0])
    _model.rotate(gamma=-90)
    return _model
# ---END-Four chamber view-----------------------------------------------------------------------------------
