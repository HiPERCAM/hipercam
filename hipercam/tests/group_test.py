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
        self.val1 = 1.
        wind = Windat(win, self.val1*np.ones((win.ny,win.nx)))
        self.agroup[1] = wind
        win = Window(100,200,10,10,1,1)
        self.val2 = 10.
        wind = Windat(win, self.val2*np.ones((win.ny,win.nx)))
        self.agroup[2] = wind

    def test_group_copy(self):
        agroup = self.agroup.copy()
        agroup[1].data[0,0] += 5
        self.assertEqual(agroup[1].data[0,0], self.agroup[1].data[0,0]+5.,
                         'Agroup copy failed')

    def test_group_setitem(self):
        group = Group()
        window = Window(1,1,100,100,1,1)
        self.assertRaises(KeyError, group.__setitem__, 'a', window)
        group[1] = window
        self.assertRaises(ValueError, group.__setitem__, 2, window)

    def test_agroup_iadd(self):
        self.agroup += 9.
        self.assertTrue(
            self.agroup[1].data[0,0] == self.val1 + 9. and
            self.agroup[2].data[0,0] == self.val2 + 9.,
            'in place addition to an Agroup has failed')

    def test_agroup_isub(self):
        self.agroup -= 9.
        self.assertTrue(
            self.agroup[1].data[0,0] == self.val1 - 9. and
            self.agroup[2].data[0,0] == self.val2 - 9.,
            'in place subtraction from an Agroup has failed')

    def test_agroup_add(self):
        agroup = self.agroup + 9.
        self.assertTrue(
            agroup[1].data[0,0] == self.val1 + 9. and
            agroup[2].data[0,0] == self.val2 + 9.,
            'addition to an Agroup has failed')

if __name__ == '__main__':
    unittest.main()
