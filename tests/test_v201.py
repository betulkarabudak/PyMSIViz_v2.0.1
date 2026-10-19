import unittest
import numpy as np
from core.volcano_analytics import select_volcano_features, build_global_spectrum_for_features

class FakeParser:
    coordinates = [(i+1,1,1) for i in range(30)]
    def getspectrum(self, i):
        return np.arange(11, dtype=float), np.ones(11) * (i + 1)

class TestVolcanoV202(unittest.TestCase):
    def test_all_zero(self):
        self.assertEqual(len(select_volcano_features(np.arange(100.), np.zeros(100))), 0)
    def test_flat_positive(self):
        self.assertEqual(len(select_volcano_features(np.arange(100.), np.full(100, 8.))), 0)
    def test_real_peaks_ordered_by_height(self):
        mz = np.arange(101., dtype=float)
        vals = np.zeros(101); vals[25] = 50; vals[65] = 100
        feats = select_volcano_features(mz, vals, max_features=2)
        self.assertGreaterEqual(len(feats), 1)
        self.assertEqual(feats[0], 65)
    def test_global_spectrum_seed(self):
        p = FakeParser()
        a = build_global_spectrum_for_features(p, max_pixels=10, random_state=9)
        b = build_global_spectrum_for_features(p, max_pixels=10, random_state=9)
        np.testing.assert_array_equal(a[1], b[1])

if __name__ == '__main__': unittest.main()

class TestROIGuards(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sys, types
        if 'pyimzml' not in sys.modules:
            pkg=types.ModuleType('pyimzml'); sub=types.ModuleType('pyimzml.ImzMLParser')
            sub.ImzMLParser=object
            sys.modules['pyimzml']=pkg;sys.modules['pyimzml.ImzMLParser']=sub
        from core.msi_tasks import task_roi_volcano
        cls.run_task=staticmethod(task_roi_volcano)
    def test_mismatched_plane(self):
        class P:
            imzmldict={'max count of pixels x':2,'max count of pixels y':2}
        mask=np.ones((2,2),dtype=bool)
        rois={'a':{'mask':mask,'z':1},'b':{'mask':mask,'z':2}}
        out=self.run_task(P(),rois,'a','b',.1)
        self.assertEqual(out['reason'],'roi_z_mismatch')
    def test_mismatched_shape(self):
        class P:
            imzmldict={'max count of pixels x':2,'max count of pixels y':2}
        rois={'a':{'mask':np.ones((2,2),dtype=bool),'z':1},'b':{'mask':np.ones((3,2),dtype=bool),'z':1}}
        out=self.run_task(P(),rois,'a','b',.1)
        self.assertEqual(out['reason'],'roi_mask_shape_mismatch')
