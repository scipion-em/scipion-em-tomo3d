# **************************************************************************
# *
# * Authors:     Scipion Team
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
from enum import Enum
from os.path import join

import mrcfile
import numpy as np
from pyworkflow.utils import makePath

from tomo.protocols import ProtTomoBase

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import IntParam, EnumParam, PointerParam, BooleanParam

from tomo.objects import Tomogram, SetOfTomograms


class outputTomoRecObjects(Enum):
    tomograms = SetOfTomograms


class ProtBaseReconstruct(EMProtocol, ProtTomoBase):
    """ Reconstruct tomograms from aligned tilt series using TOMO3D from
    Software from: https://sites.google.com/site/3demimageprocessing/
    Returns the set of tomograms
    """

    # --------------------------- DEFINE param functions --------------------------------------------
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.outputFiles = []

    def _defineParams(self, form):
        pass

    @staticmethod
    def _defineInputParams(form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label='Input')
        form.addParam('inputSetOfTiltSeries', PointerParam, important=True,
                      pointerClass='SetOfTiltSeries',
                      label='Tilt Series')

    @staticmethod
    def _defineSetShapeParams(form):
        form.addParam('setShape', BooleanParam,
                      default=True,
                      label='Set manual tomogram shape',
                      display=EnumParam.DISPLAY_HLIST,
                      help='By deafault the shape of the tomogram is defined by the tilt series shape')

        group = form.addGroup('Tomogram shape', condition='setShape')
        group.addParam('width', IntParam,
                       default=0,
                       label='Width',
                       help='Focus the tomogram in a region of the tilt series')
        group.addParam('height', IntParam,
                       default=0,
                       label='Thickness',
                       help='Height of the reconstructed tomogram (Default: width of the tomogram)')
        group.addParam('iniSlice', IntParam,
                       default=0,
                       label='Initial slice',
                       help='Initial slice (of range) to include')
        group.addParam('finSlice', IntParam,
                       default=0,
                       label='Final slice',
                       help='Final slice (of range) to include (Maximum must be the size of tilt series)')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        pass

    def rotXTomo(self, tsId):
        """Result of the reconstruction must be rotated 90 degrees around the X
        axis to recover the original orientation (due to jjsoft design)"""
        outPath = self._getExtraPath(tsId)
        makePath(outPath)
        inTomoFile = join(self._getTmpPath(tsId), '%s.mrc' % tsId)
        outTomoFile = join(outPath, '%s.mrc' % tsId)

        with mrcfile.open(inTomoFile, mode='r', permissive=True) as mrc:
            rotData = np.rot90(mrc.data)

        with mrcfile.open(outTomoFile, mode='w+') as mrc:
            mrc.set_data(rotData)

        return outTomoFile

    def createOutputStep(self):
        outputTomos = self._createSetOfTomograms()
        outputTomos.copyInfo(self.inputSetOfTiltSeries.get())

        for tomoPath, ts in zip(self.outputFiles, self.inputSetOfTiltSeries.get()):
            tomo = Tomogram()
            tomo.setLocation(tomoPath)
            tomo.setSamplingRate(ts.getSamplingRate())
            tomo.setOrigin()
            tomo.setTsId(ts.getTsId())
            outputTomos.append(tomo)

        self._defineOutputs(**{outputTomoRecObjects.tomograms.name:outputTomos})
        self._defineSourceRelation(self.inputSetOfTiltSeries, outputTomos)

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        pass

    def _methods(self):
        pass

    def _citations(self):
        return ['Fernandez2018','Fernandez2009']

    # --------------------------- UTILS functions --------------------------------------------
    @staticmethod
    def getTsFiles(tsFolder, tsId):
        """Returns the path of the Tilt Serie and the angles files"""
        prefix = join(tsFolder, tsId)
        TsPath = prefix + '.st'
        AnglesPath = prefix + '.rawtlt'
        return TsPath, AnglesPath

    def getWorkingDirName(self, tsId):
        return self._getTmpPath(tsId)


