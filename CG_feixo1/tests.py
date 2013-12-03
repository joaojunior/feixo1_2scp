import unittest

from main import transform_incidence2column

class TestTransformIncidence2Column(unittest.TestCase):
    def test_transform_incidence2column(self):
        incidence_matrix = [[1,0,1,0,0,1], [0,1,0,1,1,0], [1,1,1,0,0,0], [0,0,1,0,1,0]]
        expect = [[1,0,1,0], [0,1,1,0],[1,0,1,1],[0,1,0,0],[0,1,0,1],[1,0,0,0]]
        self.assertListEqual(expect, transform_incidence2column(incidence_matrix))
        
if __name__ == '__main__':
    unittest.main()