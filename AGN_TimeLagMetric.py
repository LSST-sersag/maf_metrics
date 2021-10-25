import os
import numpy as np
import pandas as pd
import healpy as hp
import math 
import time
from matplotlib import cm

# lsst libraries
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.db as db
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as mb


class AGN_TimeLagMetric:
    bundle = None
    lag = None
    z = None
    opsim = None
    cmap = None
    
    def __getOpsimData(self, opsim, band, name, nside = 32, outDir = 'TmpDir'):
        
        opsdb = db.OpsimDatabase(opsim)
        resultsDb = db.ResultsDb(outDir=outDir)
        metric=metrics.PassMetric(cols=['observationStartMJD', 'filter'])
        slicer = slicers.HealpixSlicer(nside)
        sqlconstraint = 'filter = \'' + band + '\''
        bundle = mb.MetricBundle(metric, slicer, sqlconstraint, runName=name)
        bgroup = mb.MetricBundleGroup({0: bundle}, opsdb, outDir=outDir, resultsDb=resultsDb)
        bgroup.runAll()
        return bundle
    
    
    def __getNquistValue(self, caden, lag, z):
        return lag/((1+z)*caden)
    
    
    
    def __init__(self, opsim, name, band, lag, z = 1, nside = 32):
        self.lag = lag
        self.opsim = opsim
        self.z = z
        self.name = name
        self.bundle = self.__getOpsimData(opsim, band, name, nside)
   
    def runAll(self):
        nquist = self.__getData(self.bundle, self.lag, self.z)
        self.__getPlots(nquist, self.name)
    
    
    def __getData(self, bundle, lag, z = 1):
        result = { 'nquist': [] }
        n = len(bundle.metricValues)
        data = bundle.metricValues.filled(0)
        for i in range(n):
            if data[i] == 0:
                result['nquist'].append(np.nan)
                continue
            mv = bundle.metricValues[i]['observationStartMJD']
            mv = np.sort(mv)
            val = (np.mean(np.diff(mv)))


            if  math.isnan(val) :
                result['nquist'].append(np.nan)
            else:
                result['nquist'].append(self.__getNquistValue(val, lag, z))

        return result
    
    def setCmap(self, cmap):
        self.cmap = cmap
    
    def __getPlots(self,data, opsim, maxval = 30):
        nquist = np.abs(np.array(data['nquist']))
        nquist[nquist<2.2] = np.nan
        hp.mollview(
            nquist,
            title=  opsim,
            unit='Threshold', cmap=self.cmap)
