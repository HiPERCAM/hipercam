import unittest
import copy

import numpy as np
from hipercam import Header, HipercamError

class TestHeader(unittest.TestCase):
    """
    Provides unit tests of the Header class
    """

    def setUp(self):
        self.head = Header()
        self.head['abc'] = (1.23, 'a comment')
        self.head['NREC'] = (1234, 'b comment')
        self.head['NRC'] = (12345, 'b comment')

    def test_getitem(self):
        self.assertEqual(
            self.head['NREC'], 1234,
            'could not recover stored header value'
        )

    def test_del(self):
        del self.head['NREC']
        self.assertEqual(
            self.head['NRC'], 12345,
            'could not recover stored header value after del'
        )

if __name__ == '__main__':
    unittest.main()
