#!/usr/bin/env python
import tempfile

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
from os.path import exists, join

from imod.constants import OUTPUT_TILTSERIES_NAME
from imod.protocols import ProtImodImportTransformationMatrix, ProtImodTsNormalization, ProtImodExcludeViews
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer
from tomo3d.protocols.protocol_base_reconstruct import outputTomoRecObjects
from tomo3d.protocols.protocol_reconstruct_tomogram import ProtTomo3dReconstrucTomo, WBP, SIRT
from tomo.protocols.protocol_ts_import import ProtImportTs

TS_03 = 'TS_03'
TS_54 = 'TS_54'


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
        print("\n", magentaStr(" Test tomo3D reconstruction with WBP".center(75, '-')))

        # preparing and launching the protocol
        ptomo3D = self.newProtocol(ProtTomo3dReconstrucTomo,
                                   inputSetOfTiltSeries=self.setOfTs,
                                   method=0)
        self.launchProtocol(ptomo3D, wait=True)
        setOfReconstructedTomograms = getattr(ptomo3D, outputTomoRecObjects.tomograms.name, None)

        # check results
        self._checkResults(setOfReconstructedTomograms)

    def testReconstructionSIRT(self):
        print("\n", magentaStr(" Test tomo3D reconstruction with SIRT".center(75, '-')))

        # preparing and launching the protocol
        ptomo3D = self.newProtocol(ProtTomo3dReconstrucTomo,
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


class TestTomogramRec(TestBaseCentralizedLayer):
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value
    nTomos = 2
    tomoWidth = 300
    bin4 = 4
    bin4SRate = DataSetRe4STATuto.unbinnedPixSize.value * 4
    tomoDimsThk300 = [960, 928, 300]

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()
        cls.tsWithAlignment = cls._runImportTrMatrix()
        cls.tsWithAliBin4 = cls._runTsPreprocess()

    @classmethod
    def _runImportTs(cls):
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
                                       anglesFrom=2,  # From tlt file
                                       voltage=DataSetRe4STATuto.voltage.value,
                                       magnification=DataSetRe4STATuto.magnification.value,
                                       sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                       amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                       samplingRate=cls.unbinnedSRate,
                                       doseInitial=DataSetRe4STATuto.initialDose.value,
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protImportTs)
        tsImported = getattr(protImportTs, 'outputTiltSeries', None)
        return tsImported

    @classmethod
    def _runImportTrMatrix(cls):
        print(magentaStr("\n==> Importing the TS' transformation matrices with IMOD:"))
        protImportTrMatrix = cls.newProtocol(ProtImodImportTransformationMatrix,
                                             filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                             filesPattern=DataSetRe4STATuto.transformPattern.value,
                                             inputSetOfTiltSeries=cls.importedTs)
        cls.launchProtocol(protImportTrMatrix)
        outTsSet = getattr(protImportTrMatrix, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runTsPreprocess(cls):
        print(magentaStr(f"\n==> Binning the TS with using a binning factor of {cls.bin4}:"))
        protTsPreprocess = cls.newProtocol(ProtImodTsNormalization,
                                           inputSetOfTiltSeries=cls.tsWithAlignment,
                                           binning=cls.bin4)
        cls.launchProtocol(protTsPreprocess)
        outTsSet = getattr(protTsPreprocess, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runExcludeViews(cls, inTsSet=None):
        print(magentaStr("\n==> Excluding some views with IMOD:"))
        excViewsFileName = join('/tmp', 'excViews.txt')
        with open(excViewsFileName, mode='w') as eFile:
            eFile.write(f"{TS_03} 1-3,39,40\n"
                        f"{TS_54} 1,2,39-41")
        excViewsProt = cls.newProtocol(ProtImodExcludeViews,
                                       inputSetOfTiltSeries=inTsSet,
                                       excludeViewsFile=excViewsFileName)
        cls.launchProtocol(excViewsProt)
        outTsSet = getattr(excViewsProt, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runTomoRec(cls, inTsSet=None, recMethod=None, excludedViews=False):
        method = 'WBP' if recMethod is WBP else 'SIRT'
        excViewsMsg = 'with excViews' if excludedViews else ''
        print(magentaStr(f"\n==> Reconstructing the tomograms using the method {method} {excViewsMsg}:"))
        protTomoRec = cls.newProtocol(ProtTomo3dReconstrucTomo,
                                      inputSetOfTiltSeries=inTsSet,
                                      method=recMethod,
                                      height=cls.tomoWidth,
                                      nIterations=10,
                                      numberOfThreads=8)
        protTomoRec.setObjLabel(f'Tomo rec {method} {excViewsMsg}')
        cls.launchProtocol(protTomoRec)
        outTomos = getattr(protTomoRec, outputTomoRecObjects.tomograms.name, None)
        return outTomos

    def testRecWbp(self):
        recTomos = self._runTomoRec(recMethod=WBP, inTsSet=self.tsWithAliBin4)
        self._checkTomos(recTomos)

    def testRecWbpExcludedViews(self):
        tsSetWithExcViews = self._runExcludeViews(inTsSet=self.tsWithAliBin4)
        recTomos = self._runTomoRec(recMethod=WBP, inTsSet=tsSetWithExcViews, excludedViews=True)
        self._checkTomos(recTomos)

    def testRecSirt(self):
        recTomos = self._runTomoRec(recMethod=SIRT, inTsSet=self.tsWithAliBin4)
        self._checkTomos(recTomos)

    def _checkTomos(self, inTomoSet):
        self.checkTomograms(inTomoSet,
                            expectedSetSize=self.nTomos,
                            expectedSRate=self.bin4SRate,
                            expectedDimensions=self.tomoDimsThk300)
