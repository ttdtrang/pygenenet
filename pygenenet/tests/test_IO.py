from LearningCausalNetwork import DataIO, Experiment
import numpy as np
import unittest

class DataIOTest(unittest.TestCase):
    tsd_file = 'test_data/prob1/run-1.tsd'
    def test_1_read_tsd(self):
        h, d = DataIO.read_TimeSeries(self.tsd_file, format='sequential')
        self.assertEqual(len(h) , 5)
        self.assertEqual(len(d) , 211)

class ExperimentTest(unittest.TestCase):
    tsd_file = 'test_data/prob1/run-1.tsd'
    def test_1_init(self):
        tsd = Experiment.TimeSeriesData(self.tsd_file, format='sequential')
        # self.assertEqual(len(tsd.columns) , 5)
        self.assertEqual(len(tsd) , 211)
        self.assertTrue((tsd.columns == ["time","CI","PR","CII","PRE"]).all())
        self.assertTrue(( tsd['time'] == np.arange(0, 2101, 10)).all() )
        self.assertEqual( tsd.loc[20,'CII'], 16)
    def test_2_to_csv(self):
        tsd = Experiment.TimeSeriesData(self.tsd_file, format='sequential')
        tsd.to_csv('test_output.csv') 
    def test_3_snapshot(self):
        tsd = Experiment.TimeSeriesData(self.tsd_file, format='sequential')
        t50 = tsd.snapshot(50)
        self.assertEqual(t50['CI'], 10)
        self.assertEqual(t50['PR'], 2)
        self.assertEqual(t50['CII'], 24)
        self.assertEqual(t50['PRE'], 2)
    def test_4_expset_fromfiles(self):
        exp_set = Experiment.ExperimentSet(['test_data/prob1/run-1.tsd','test_data/prob1/run-2.tsd' ], format='sequential')
        self.assertTrue((exp_set.datapoint(0, 50, 'CII') ==  24).all() )
        self.assertTrue((exp_set.datapoint(0, 50, 'CI') ==  10).all() )
        self.assertTrue((exp_set.species() ==  ['CI',  'PR', 'CII', 'PRE']) )
    def test_5_expset_fromobjs(self):
        pass
        # tsd1 = Experiment.TimeSeriesData('test_data/prob1/run-1.tsd', format='sequential')
        # tsd2 = Experiment.TimeSeriesData('test_data/prob1/run-2.tsd', format='sequential')
        # exp_set = Experiment.ExperimentSet([tsd1, tsd2 ])
    # def test_3_to_binary(self):
    #     tsd = Experiment.TimeSeriesData(self.tsd_file, format='sequential')
    #     tsd.to_binary('test_output.npy') 
    
def suite():
    suite1 = unittest.TestLoader().loadTestsFromTestCase(DataIOTest)
    suite2 =  unittest.TestLoader().loadTestsFromTestCase(ExperimentTest)
    #suite = unittest.TestSuite([suite1, suite2])
    suite = unittest.TestSuite([suite1, suite2])
    return suite

if __name__ == '__main__':
    import sys
    sys.path.append('/media/tk/DATA/Users/Trang/courses_u/ModelingBioNetworks/Project/')
    unittest.main()

    exp_set = Experiment.ExperimentSet(['test_data/prob1/run-1.tsd','test_data/prob1/run-2.tsd' ], format='sequential')
    bins = np.linspace(1,100,1)
    inds = np.digitize(exp_set['CI'], bins)
    # print exp_set[('experiment' == 1) & ('time' == 100) & ('species' == 'PRE') ]
