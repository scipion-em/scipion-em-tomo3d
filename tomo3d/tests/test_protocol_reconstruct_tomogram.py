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
import numpy as np
from os.path import exists
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo3d.protocols.protocol_base_reconstruct import outputTomoRecObjects
from tomo3d.protocols.protocol_reconstruct_tomogram import ProtJjsoftReconstructTomogram
from tomo.protocols.protocol_ts_import import ProtImportTs


class TestTomogramReconstruction(BaseTest):

    jjsoftDataTest = None
    getFile = None
    samplingRate = 1.35

    @classmethod
    def setUpClass(cls):
        """Prepare the data that we will use later on."""
        setupTestProject(cls)  # defined in BaseTest, creates cls.proj
        cls.jjsoftDataTest = DataSet.getDataSet('tomo-em')
        cls.getFile = cls.jjsoftDataTest.getFile('etomo')
        cls.setOfTs = cls._runImportTiltSeries()

    @classmethod
    def _runImportTiltSeries(cls):
        protImport = cls.newProtocol(
            ProtImportTs,
            filesPath=cls.getFile,
            filesPattern='BB{TS}.st',
            minAngle=-55,
            maxAngle=65,
            stepAngle=2,
            voltage=300,
            magnification=105000,
            sphericalAberration=2.7,
            amplitudeContrast=0.1,
            samplingRate=cls.samplingRate,
            doseInitial=0,
            dosePerFrame=0.3,
            tiltAxisAngle=90)
        cls.launchProtocol(protImport, wait=True)
        return protImport.outputTiltSeries

    def testReconstructionWBP(self):
        print ("\n", magentaStr(" Test tomo3D reconstruction with WBP".center(75, '-')))

        # preparing and launching the protocol
        ptomo3D = self.newProtocol(ProtJjsoftReconstructTomogram,
                                   inputSetOfTiltSeries=self.setOfTs,
                                   method=0)
        self.launchProtocol(ptomo3D, wait=True)
        setOfReconstructedTomograms = getattr(ptomo3D, outputTomoRecObjects.tomograms.name, None)

        # check results
        self._checkResults(setOfReconstructedTomograms)

    def testReconstructionSIRT(self):
        print ("\n", magentaStr(" Test tomo3D reconstruction with SIRT".center(75, '-')))

        # preparing and launching the protocol
        ptomo3D = self.newProtocol(ProtJjsoftReconstructTomogram,
                                   inputSetOfTiltSeries=self.setOfTs,
                                   method=1,
                                   nIterations=1)
        self.launchProtocol(ptomo3D, wait=True)
        setOfReconstructedTomograms = getattr(ptomo3D, outputTomoRecObjects.tomograms.name, None)

        # check results
        self._checkResults(setOfReconstructedTomograms)

    def _checkResults(self, outputSet):
        # Set of tomograms checks
        self.assertIsNotNone(outputSet, "There was some problem with the output")
        self.assertEqual(outputSet.getSize(), self.setOfTs.getSize(), "The number of the denoised tomograms is wrong")
        self.assertEqual(outputSet.getSamplingRate(), self.samplingRate)
        # Tomograms checks
        testTsIds = ['b', 'a']
        testOrigin = np.array([
            [1.0, 0.0, 0.0, -345.6],
            [0.0, 1.0, 0.0, -345.6],
            [0.0, 0.0, 1.0, -345.6],
            [0.0, 0.0, 0.0, 1.0]])
        for i, tomo in enumerate(outputSet):
            self.assertTrue(exists(tomo.getFileName()))
            self.assertEqual(tomo.getTsId(), testTsIds[i])
            self.assertEqual(tomo.getSamplingRate(), self.samplingRate)
            self.assertTrue(np.allclose(tomo.getOrigin(force=True).getMatrix(), testOrigin, rtol=1e-2))


