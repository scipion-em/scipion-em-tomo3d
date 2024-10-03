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
from pyworkflow.object import Set
from tomo.protocols import ProtTomoBase
from pwem.protocols import EMProtocol
from tomo.objects import Tomogram, SetOfTomograms

# Odd/even
EVEN = 'even'
ODD = 'odd'


class outputTomo3dObjects(Enum):
    tomograms = SetOfTomograms


class ProtBaseTomo3d(EMProtocol, ProtTomoBase):
    """ Reconstruct tomograms from aligned tilt series using TOMO3D from
    Software from: https://sites.google.com/site/3demimageprocessing/
    Returns the set of tomograms
    """
    _OUTNAME = outputTomo3dObjects.tomograms.name
    _possibleOutputs = {_OUTNAME: SetOfTomograms}

    # --------------------------- DEFINE param functions --------------------------------------------
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.objDict = None  # {tsId: obj}, where obj can be a TS or a tomogram

    def _defineParams(self, form):
        pass

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        pass

    def rotXTomo(self, tsId, suffix=None):
        """Result of the reconstruction must be rotated 90 degrees around the X
        axis to recover the original orientation (due to tomo3d design)"""

        inTomoFile = self._getTmpTomoOutFName(tsId, suffix=suffix)
        outTomoFile = self._getOutTomoFile(tsId, suffix=suffix)

        with mrcfile.mmap(inTomoFile, mode='r', permissive=True) as mrc:
            rotData = np.rot90(mrc.data)

        with mrcfile.mmap(outTomoFile, mode='w+') as mrc:
            mrc.set_data(rotData)

    def createOutputStep(self, tsId, doEvenOdd=False):
        obj = self.objDict[tsId]
        acq = obj.getAcquisition().clone()
        outputTomos = self._getOutputSetOfTomograms()
        tomo = Tomogram()
        tomo.setTsId(obj.getTsId())
        tomo.setLocation(self._getOutTomoFile(tsId))
        tomo.setSamplingRate(obj.getSamplingRate())
        tomo.setAcquisition(acq)
        tomo.setOrigin()
        if doEvenOdd:
            tomo.setHalfMaps([self._getOutTomoFile(tsId, suffix=EVEN), self._getOutTomoFile(tsId, suffix=ODD)])
        outputTomos.append(tomo)
        outputTomos.update(tomo)
        outputTomos.write()
        self._store(outputTomos)

    def closeOutputSetsStep(self):
        self._closeOutputSet()

    # --------------------------- INFO functions --------------------------------------------

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

    def _getOutputSetOfTomograms(self):
        outTomograms = getattr(self, outputTomo3dObjects.tomograms.name, None)
        if outTomograms:
            outTomograms.enableAppend()
            tomograms = outTomograms
        else:
            inSet = self.getInputSet()
            tomograms = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            tomograms.copyInfo(inSet)
            tomograms.setStreamState(Set.STREAM_OPEN)
            setattr(self, outputTomo3dObjects.tomograms.name, tomograms)
            self._defineOutputs(**{self._OUTNAME: tomograms})

            self._defineSourceRelation(inSet, tomograms)

        return tomograms

    def _getOutTomoFile(self, tsId, suffix=None):
        return join(self._getTsExtraDir(tsId), f'{tsId}{self._manageSuffix(suffix)}.mrc')

    def _getTmpTomoOutFName(self, tsId, suffix=None):
        return join(self._getTsTmpDir(tsId), f'{tsId}{self._manageSuffix(suffix)}.mrc')

    @staticmethod
    def _manageSuffix(suffix):
        return f'_{suffix}' if suffix else ''

    def _getTsTmpDir(self, tsId):
        return self._getTmpPath(tsId)

    def _getTsExtraDir(self, tsId):
        return self._getExtraPath(tsId)

    def getInputSet(self):
        pass
