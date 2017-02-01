import unittest
import copy
from astropy.io import fits

import numpy as np
from hipercam import Window, Windat, Group, Agroup, HipercamError

class TestGroup(unittest.TestCase):
    """Provides simple tests of Group methods

    """

    def setUp(self):
        self.agroup = Agroup()
        win = Window(1,1,10,10,1,1)
        wind = Windat(win, np.ones((win.ny,win.nx)))
        self.agroup[1] = wind
        win = Window(100,200,10,10,1,1)
        wind = Windat(win, 10.*np.ones((win.ny,win.nx)))
        self.agroup[2] = wind

    def test_group_setitem(self):
        group = Group()
        window = Window(1,1,100,100,1,1)
        self.assertRaises(HipercamError, group.__setitem__, 'a', window)
        group[1] = window
        self.assertRaises(HipercamError, group.__setitem__, 2, window)

    def test_agroup_iadd(self):
        agroup = copy.deepcopy(self.agroup)
        agroup += 9.
        self.assertTrue(
            agroup[1].data[0,0] == 10. and agroup[2].data[0,0] == 19.,
            'addition to an Agroup has failed')

    def test_agroup_isub(self):
        agroup = copy.deepcopy(self.agroup)
        agroup -= 9.
        self.assertTrue(
            agroup[1].data[0,0] == -8. and agroup[2].data[0,0] == 1.,
            'subtraction from an Agroup has failed')

if __name__ == '__main__':
    unittest.main()
