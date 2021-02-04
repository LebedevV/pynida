from pynida.backend import oneshot_contours
from pynida.simple_functions import get_n_frame
import cv2
import pickle
import numpy as np
import os

package_directory = os.path.dirname(os.path.abspath(__file__))

# Function used to generate test frames from video.
# Video itself is not included to the repo due to size considerations.
def pickle_frame(path_to_video, frame_n, savepath):
    frame = get_n_frame(path_to_video, 1, 100, frame_n)
    frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    with open(savepath, "wb") as file:
        pickle.dump(frame, file)


def test_oneshot_contours_frame():
    path = './data/testframe.pkl'
    frame = pickle.load(open(os.path.join(package_directory, path), "rb"))
    results = oneshot_contours(frame, None, False)
    expected_results = {
        'mid_x': 1073.0,
        'mid_y': 11.0,
        'ts_len': 1348.8087336609294,
        'invvol_t': 0.003989385581130855,
        'sp_len': 1299.0400301761297,
        'mid_len': 1073.0563824888234,
        'vol_t': 37287325.57218897,
        'vol_b': 1282516737.3145592,
        'invvol_b': 0.0032679097890607653,
        'corrvol_b': 1620777.6342996806,
        'corrvol_t': -125917908.51512668
    }
    for key, value in expected_results.items():
        assert np.isclose(results[key], value), f"{key} is different: got {results[key]}, expected {value}"
