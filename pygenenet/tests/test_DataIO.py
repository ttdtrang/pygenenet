from pygenenet import DataIO, Experiment
import numpy as np
import unittest

class DataIOTest(unittest.TestCase):
    tsd_file = 'test_data/prob1/run-1.tsd'
    def test_1_read_tsd(self):
        h, d = DataIO.read_TimeSeries(self.tsd_file, format='sequential')
        self.assertEqual(len(h) , 5)
        self.assertEqual(len(d) , 211)

class ExperimentTest(unittest.TestCase):
    def test_1_init(self):
        tsd = DataIO.TimeSeriesData(self.tsd_file, format='sequential')
        self.assertEqual(len(tsd.columns) , 5)
        self.assertEqual(len(tsd.data) , 211)
        self.assertEqual(len(tsd['time']) , 211)
        self.assertTrue(( tsd['time'] == np.arange(0, 2101, 10)).all() )
        self.assertEqual( tsd['CII'][2] , 16)
    def test_2_to_csv(self):
        tsd = DataIO.TimeSeriesData(self.tsd_file, format='sequential')
        tsd.to_csv('test_output.csv') 
    def test_3_to_binary(self):
        tsd = DataIO.TimeSeriesData(self.tsd_file, format='sequential')
        tsd.to_binary('test_output.npy') 

def suite():
    suite1 = unittest.TestLoader().loadTestsFromTestCase(DataIOTest)
    suite2 =  unittest.TestLoader().loadTestsFromTestCase(ExperimentTest)
    suite = unittest.TestSuite([suite1, suite2])
    return suite

if __name__ == '__main__':
    import sys
    sys.path.append('/media/tk/DATA/Users/Trang/courses_u/ModelingBioNetworks/Project/')
    unittest.main()
