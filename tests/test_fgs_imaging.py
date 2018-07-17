
import os

from mirage.scripts import imaging_simulator as im

os.environ['TEST_DATA'] = os.path.join(os.path.dirname(__file__), 'test_data/FGS')

def test_fgs_imaging():
    m = im.ImgSim()
    m.paramfile = os.path.join(os.path.dirname(__file__), 'test_data/FGS/fgs_imaging_example.yaml')
    m.create()

