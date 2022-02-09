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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import greenStr
from imod.protocols import ProtImodXcorrPrealignment, ProtImodFiducialAlignment

from tomo3d.protocols.protocol_align_ts import ProtJjsoftAlignTs
from tomo3d.protocols.protocol_motion_compensation import ProtJjsoftAlignReconstructTomogram

from tomo.protocols.protocol_ts_import import ProtImportTs


class TestAlignTs(BaseTest):

    @classmethod
    def setUpClass(cls):
        """Prepare the data that we will use later on."""
        setupTestProject(cls)  # defined in BaseTest, creates cls.proj
        cls.jjsoftDataTest = DataSet.getDataSet('tomo-em')
        cls.getFile = cls.jjsoftDataTest.getFile('etomo')
        cls.setOfTs = cls._runImportTiltSeries()
        cls.setOfXcorrTs = cls._runXcorrPrealignment()
        cls.setOfFiducials, cls.setOfFiducialTs = cls._runFiducialModel()

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
            samplingRate=1.35,
            doseInitial=0,
            dosePerFrame=0.3,
            tiltAxisAngle=90
        )
        cls.launchProtocol(protImport, wait=True)
        return protImport.outputTiltSeries

    @classmethod
    def _runXcorrPrealignment(cls):
        protXcorr = cls.newProtocol(ProtImodXcorrPrealignment,
                                    inputSetOfTiltSeries=cls.setOfTs,
                                    computeAlignment=1,
                                    binning=1,
                                    rotationAngle=0.0)
        cls.launchProtocol(protXcorr)
        return protXcorr.outputSetOfTiltSeries

    @classmethod
    def _runFiducialModel(cls):
        protFiducials = cls.newProtocol(ProtImodFiducialAlignment,
                                        inputSetOfTiltSeries=cls.setOfXcorrTs,
                                        fiducialDiameter=4.95,
                                        binning=1)
        cls.launchProtocol(protFiducials)
        return protFiducials.outputFiducialModelNoGaps, protFiducials.outputSetOfTiltSeries

    # The tests themselves.
    #
    def testWarpAlign(self):
        print ("\n", greenStr(" Test Ts alignment with warpalign ".center(75, '-')))

        # preparing and launching the protocol
        pwarpalign = self.newProtocol(ProtJjsoftAlignTs,
                                      inputSetOfTiltSeries=self.setOfFiducialTs,
                                      inputSetOfLandmarkModels=self.setOfFiducials)
        self.launchProtocol(pwarpalign, wait=True)
        setOfAlignedTs = pwarpalign.outputInterpolatedSetOfTiltSeries

        # some general assertions
        self.assertIsNotNone(setOfAlignedTs,
                             "There was some problem with the output")
        self.assertEqual(setOfAlignedTs.getSize(), self.setOfFiducialTs.getSize(),
                         "The number of the aligned Ts is wrong")

    def testAlignAndReconstruct(self):
        print("\n", greenStr(" Test Ts alignment and tomogram reconstruction with Jjsof".center(75, '-')))

        # preparing and launching the protocol
        palireco = self.newProtocol(ProtJjsoftAlignReconstructTomogram,
                                    inputSetOfTiltSeries=self.setOfFiducialTs,
                                    inputSetOfLandmarkModels=self.setOfFiducials)
        self.launchProtocol(palireco, wait=True)
        setOfTomograms = palireco.outputTomograms

        # some general assertions
        self.assertIsNotNone(setOfTomograms,
                             "There was some problem with the output")
        self.assertEqual(setOfTomograms.getSize(), self.setOfFiducialTs.getSize(),
                         "The number of tomograms is wrong")
