#!/usr/bin/env python

# ***************************************************************************
# *
# * Authors:     Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from os.path import exists

import numpy as np
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr, removeBaseExt
from tomo3d.protocols.protocol_denoise_tomogram import ProtJjsoftProtDenoiseTomogram, DENOISE_EED, DENOISE_BF, \
    outputDenoiseObjects
from tomo.protocols.protocol_import_tomograms import ProtImportTomograms


class TestTomogramDenoising(BaseTest):

    jjsoftDataTest = None
    setOfTomograms = None
    samplingRate = 1.35

    @classmethod
    def setUpClass(cls):
        """Prepare the data that we will use later on."""
        setupTestProject(cls)  # defined in BaseTest, creates cls.proj

        cls.jjsoftDataTest = DataSet.getDataSet('tomo-em')
        cls.setOfTomograms = cls._importTomograms()

    @classmethod
    def _importTomograms(cls):
        """ Importing a set of tomograms
        """
        pImpTomograms = cls.newProtocol(ProtImportTomograms,
                                        filesPath=cls.jjsoftDataTest.getFile('tomo'),
                                        samplingRate=cls.samplingRate,
                                        acquisitionAngleMax=40.0,
                                        acquisitionAngleMin=-40.0,
                                        tiltAxisAngle=90)

        pImpTomograms.setObjLabel('Import tomograms')

        # we launch the protocol to obtain the tomograms
        cls.launchProtocol(pImpTomograms, wait=True)

        # Setting the set of tomograms object
        return pImpTomograms.outputTomograms

    def _runDenoising(self, denoisingMethod):
        # preparing and launching the protocol
        pDenoiseEED = self.newProtocol(ProtJjsoftProtDenoiseTomogram,
                                       inputSetTomograms=self.setOfTomograms,
                                       method=denoisingMethod,
                                       SigmaGaussian=0.5,
                                       nIter=1,
                                       TimeStep=0.1,
                                       Lambda=-1.0)
        self.launchProtocol(pDenoiseEED, wait=True)
        return getattr(pDenoiseEED, outputDenoiseObjects.tomograms.name, None)

    def testDenoisingEED(self):
        print("\n", magentaStr(" Test EED denoising ".center(75, '-')))
        # preparing and launching the protocol
        setOfEEDDenoisedTomograms = self._runDenoising(DENOISE_EED)
        # check results
        self._checkResults(setOfEEDDenoisedTomograms)

    def testDenoisingBFlow(self):
        print ("\n", magentaStr(" Test BFlow denoising ".center(75, '-')))
        # preparing and launching the protocol
        setOfEEDDenoisedTomograms = self._runDenoising(DENOISE_BF)
        # check results
        self._checkResults(setOfEEDDenoisedTomograms)

    def _checkResults(self, outputSet):
        self.assertIsNotNone(outputSet, "There was some problem with the output")
        self.assertEqual(outputSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(outputSet.getSize(), self.setOfTomograms.getSize(),
                         "The number of the denoised tomograms is wrong")

        # Tomograms checks
        testOrigin = np.array([
            [1.0, 0.0, 0.0, -691.2],
            [0.0, 1.0, 0.0, -691.2],
            [0.0, 0.0, 1.0, -345.6],
            [0.0, 0.0, 0.0, 1.0]])

        for tomo in outputSet:
            self.assertTrue(exists(tomo.getFileName()))
            self.assertEqual(tomo.getTsId(), removeBaseExt(tomo.getFileName()))  # Tomograms were imported, so the tsId would be the tomo basename
            self.assertEqual(tomo.getSamplingRate(), self.samplingRate)
            self.assertTrue(np.allclose(tomo.getOrigin(force=True).getMatrix(), testOrigin, rtol=1e-2))
