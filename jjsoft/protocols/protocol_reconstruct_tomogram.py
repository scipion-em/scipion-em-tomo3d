# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (daniel.delhoyo.gomez@alumnos.upm.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from os.path import join

import mrcfile
import numpy as np

from jjsoft.protocols.protocol_base_reconstruct import ProtBaseReconstruct
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.utils import makePath

from jjsoft import Plugin
from tomo.protocols import ProtTomoBase

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import IntParam, EnumParam, PointerParam, FloatParam

from tomo.objects import Tomogram
import os


class ProtJjsoftReconstructTomogram(ProtBaseReconstruct):
    """ Reconstruct tomograms from aligned tilt series using TOMO3D from
    Software from: https://sites.google.com/site/3demimageprocessing/
    Returns the set of tomograms
    """
    _label = 'reconstruct tomogram'
    _devStatus = BETA

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineInputParams(form)
        self._defineReconstructParams(form)
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert every step of the protocol"""
        self._insertFunctionStep(self.convertInputStep)
        self.outputFiles = []
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getTmpPath(tsId)
            self._insertFunctionStep(self.reconstructTomogramStep, tsId, workingFolder)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getTmpPath(tsId)
            prefix = os.path.join(workingFolder, tsId)
            makePath(workingFolder)

            outputStackFn = prefix + '.st'
            outputTltFn = prefix + '.rawtlt'

            ts.applyTransform(outputStackFn)
            ts.generateTltFile(outputTltFn)

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        pass

    def _methods(self):
        pass

    def _citations(self):
        return ['Fernandez2018', 'Fernandez2009']

    # # --------------------------- UTILS functions --------------------------------------------
    # def get_Ts_files(self,ts_folder,TsId):
    #     """Returns the path of the Tilt Serie and the angles files"""
    #     prefix = os.path.join(ts_folder, TsId)
    #     TsPath = prefix + '.st'
    #     AnglesPath = prefix + '.rawtlt'
    #     return TsPath, AnglesPath



