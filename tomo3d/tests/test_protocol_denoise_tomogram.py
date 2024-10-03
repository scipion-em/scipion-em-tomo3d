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
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo.tests import DataSetRe4STATuto, RE4_STA_TUTO
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer
from tomo3d.protocols.protocol_base import outputTomo3dObjects
from tomo3d.protocols.protocol_denoise_tomogram import ProtTomo3dProtDenoiseTomogram, DENOISE_EED, \
    DENOISE_BF
from tomo.protocols.protocol_import_tomograms import ProtImportTomograms, OUTPUT_NAME


class TestTomoDenoising(TestBaseCentralizedLayer):
    bin4SRate = DataSetRe4STATuto.unbinnedPixSize.value * 4

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls.importedTomos = cls._runImportTomograms()
        # 5 TS with no. tilt-images:
        #   - TS_01 = 40
        #   - TS_03 = 40
        #   - TS_43 = 41
        #   - TS_45 = 41
        #   - TS_54 = 41
        #
        # 5 tomograms with a thickness of (px):
        #   - TS_01 = 340
        #   - TS_03 = 280
        #   - TS_43 = 300
        #   - TS_45 = 300
        #   - TS_54 = 280

    @classmethod
    def _runImportTomograms(cls):
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomos = cls.newProtocol(ProtImportTomograms,
                                          filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                          filesPattern='TS_0*.mrc',  # TS_01 and TS_03
                                          samplingRate=DataSetRe4STATuto.sRateBin4.value)  # Bin 4
        cls.launchProtocol(protImportTomos)
        outTomos = getattr(protImportTomos, OUTPUT_NAME, None)
        return outTomos

    @classmethod
    def _runDenoising(cls, denosingMethod=None):
        denoisingStr = 'BFlow' if denosingMethod == DENOISE_BF else 'EED'
        print(magentaStr(f"\n==> Denoising the tomograms using the method {denoisingStr}:"))
        protDenosing = cls.newProtocol(ProtTomo3dProtDenoiseTomogram,
                                       inputSetTomograms=cls.importedTomos,
                                       method=denosingMethod,
                                       numberOfThreads=8)
        cls.launchProtocol(protDenosing)
        protDenosing.setObjLabel(denoisingStr)
        return getattr(protDenosing, outputTomo3dObjects.tomograms.name, None)

    def testDenoisingBFlow(self):
        denoisedTomos = self._runDenoising(denosingMethod=DENOISE_BF)
        self._checkResults(denoisedTomos)

    def testDenoisingEED(self):
        denoisedTomos = self._runDenoising(denosingMethod=DENOISE_EED)
        self._checkResults(denoisedTomos)

    def _checkResults(self, tomoSet):
        TS_01 = 'TS_01'
        TS_03 = 'TS_03'
        tomoDimsThk280 = [928, 928, 280]
        tomoDimsThk340 = [928, 928, 340]
        expectedDimensionsDict = {
            TS_01: tomoDimsThk340,
            TS_03: tomoDimsThk280,
        }
        self.checkTomograms(tomoSet,
                            expectedSetSize=len(expectedDimensionsDict),
                            expectedSRate=self.bin4SRate,
                            expectedDimensions=expectedDimensionsDict,
                            isHeterogeneousSet=True)
