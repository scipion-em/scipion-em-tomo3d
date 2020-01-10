#!/usr/bin/env python

# ***************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
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

import random

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import redStr, greenStr

from jjsoft.protocols.protocol_denoise_tomogram import JjsoftProtDenoiseTomogram

from tomo.objects import Tomogram, SetOfTomograms
from tomo.protocols.protocol_base import ProtTomoImportFiles, ProtTomoImportAcquisition
from tomo.protocols.protocol_import_tomograms import ProtImportTomograms


class TestTomogramDenoising(BaseTest):

    @classmethod
    def setUpClass(cls):
        """Prepare the data that we will use later on."""
        setupTestProject(cls)  # defined in BaseTest, creates cls.proj

        cls.jjsoftDataTest = DataSet.getDataSet('tomo-em')

        def _importTomograms():
            """ We import a set of classes with an alterate import particles
                because there is not a import classes protocol
            """
            pImpTomograms = cls.proj.newProtocol(ProtImportTomograms,
                                               filesPath=cls.jjsoftDataTest.getFile('tomo'),
                                               samplingRate=1.35,
                                               acquisitionAngleMax=40.0,
                                               acquisitionAngleMin=-40.0)

            pImpTomograms.setObjLabel('Import tomograms')

            # we launch the protocol to obtain the tomograms
            cls.proj.launchProtocol(pImpTomograms, wait=True)

            # Setting the set of tomograms object
            tomoSet = pImpTomograms.outputTomograms

            # Defining the set of tomograms as output
            pImpTomograms._defineOutputs(outputClasses=tomoSet)

            return pImpTomograms.outputClasses

        cls.setOfTomograms = _importTomograms()


    # The tests themselves.
    #
    def testDenoisingAND(self):
        print "\n", greenStr(" Test AND denoising ".center(75, '-'))

        # preparing and launching the protocol
        pDenoiseAND = self.proj.newProtocol(JjsoftProtDenoiseTomogram,
                                            inputSetTomograms=self.setOfTomograms,
                                            method=0,
                                            SigmaGaussian=0.5,
                                            nIter=1,
                                            TimeStep=0.1)
        self.proj.launchProtocol(pDenoiseAND, wait=True)
        setOfANDDenoisedTomograms = pDenoiseAND.outputTomograms

        # some general assertions
        self.assertIsNotNone(setOfANDDenoisedTomograms,
                             "There was some problem with the output")
        self.assertEqual(setOfANDDenoisedTomograms.getSize(), self.setOfTomograms.getSize(),
                         "The number of the denoised tomograms is wrong")


    def testDenoisingEED(self):
        print "\n", greenStr(" Test EED denoising ".center(75, '-'))

        # preparing and launching the protocol
        pDenoiseEED = self.proj.newProtocol(JjsoftProtDenoiseTomogram,
                                            inputSetTomograms=self.setOfTomograms,
                                            method=2,
                                            SigmaGaussian=0.5,
                                            nIter=1,
                                            TimeStep=0.1)
        self.proj.launchProtocol(pDenoiseEED, wait=True)
        setOfEEDDenoisedTomograms = pDenoiseEED.outputTomograms

        # some general assertions
        self.assertIsNotNone(setOfEEDDenoisedTomograms,
                             "There was some problem with the output")
        self.assertEqual(setOfEEDDenoisedTomograms.getSize(), self.setOfTomograms.getSize(),
                         "The number of the denoised tomograms is wrong")

    def testDenoisingBFlow(self):
        print "\n", greenStr(" Test BFlow denoising ".center(75, '-'))

        # preparing and launching the protocol
        pDenoiseBFlow = self.proj.newProtocol(JjsoftProtDenoiseTomogram,
                                            inputSetTomograms=self.setOfTomograms,
                                            method=1,
                                            SigmaGaussian=0.5,
                                            nIter=1,
                                            TimeStep=0.1)
        self.proj.launchProtocol(pDenoiseBFlow, wait=True)
        setOfBFlowDenoisedTomograms = pDenoiseBFlow.outputTomograms

        # some general assertions
        self.assertIsNotNone(setOfBFlowDenoisedTomograms,
                             "There was some problem with the output")
        self.assertEqual(setOfBFlowDenoisedTomograms.getSize(), self.setOfTomograms.getSize(),
                         "The number of the denoised tomograms is wrong")