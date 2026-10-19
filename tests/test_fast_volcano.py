import unittest
import numpy as np
try:
    from core.msi_tasks import _volcano_roi_intensity_matrix
    from core.utils import getionimage
except ImportError:
    import sys, types
    if 'pyimzml.ImzMLParser' not in sys.modules:
        pkg=types.ModuleType('pyimzml')
        mod=types.ModuleType('pyimzml.ImzMLParser')
        mod.ImzMLParser=object
        sys.modules['pyimzml']=pkg
        sys.modules['pyimzml.ImzMLParser']=mod
    from core.msi_tasks import _volcano_roi_intensity_matrix
    from core.utils import getionimage

class FakeMSI:
    def __init__(self):
        self.imzmldict={'max count of pixels x':11,'max count of pixels y':12}
        self.coordinates=[(x,y,1) for y in range(1,13) for x in range(1,12) if (x+y)%7 != 0]
        self.reads=0
        self.mass=np.arange(100.,120.,0.1)
    def getspectrum(self,i):
        self.reads+=1
        intensities=np.sin(self.mass*3 + i*.2)**2*100 + (i%9)
        return self.mass,intensities

class TestSinglePass(unittest.TestCase):
    def test_matches_legacy_and_reduces_reads(self):
        p=FakeMSI()
        features=np.array([101.2,102.9,105.0,110.8,117.5])
        selected=np.array([0,2,5,11,24,31,44,70,88,100,131])
        tol=0.21
        new=_volcano_roi_intensity_matrix(p,features,tol,1,selected)
        new_reads=p.reads
        self.assertLessEqual(new_reads,len(selected))
        p.reads=0
        old=np.column_stack([getionimage(p,float(mz),tol,1).ravel()[selected] for mz in features])
        old_reads=p.reads
        np.testing.assert_allclose(new,old,rtol=1e-13,atol=1e-11)
        self.assertLess(new_reads,old_reads)
    def test_cancel(self):
        p=FakeMSI()
        result=_volcano_roi_intensity_matrix(p,[102.],0.2,1,[1,10],cancel_cb=lambda: True)
        self.assertIsNone(result)

if __name__=='__main__': unittest.main()
