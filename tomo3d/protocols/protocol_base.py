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
import logging
from enum import Enum
from os.path import join
from typing import Union

import mrcfile
import numpy as np
from pyworkflow.object import Set, Boolean
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import cyanStr
from tomo.protocols import ProtTomoBase
from pwem.protocols import EMProtocol
from tomo.objects import Tomogram, SetOfTomograms, SetOfTiltSeries, TiltSeries

logger = logging.getLogger(__name__)

# Odd/even
EVEN = 'even'
ODD = 'odd'
DO_EVEN_ODD = 'doEvenOdd'

# Extensions
ST_EXT = '.st'
MRC_EXT = '.mrc'
MRCS_EXT = '.mrcs'
RAWTLT_EXT = '.rawtlt'


class outputTomo3dObjects(Enum):
    tomograms = SetOfTomograms


class ProtBaseTomo3d(EMProtocol, ProtTomoBase):
    """ Reconstruct tomograms from aligned tilt series using TOMO3D from
    Software from: https://sites.google.com/site/3demimageprocessing/
    Returns the set of tomograms
    """
    _OUTNAME = outputTomo3dObjects.tomograms.name
    _possibleOutputs = {_OUTNAME: SetOfTomograms}
    stepsExecutionMode = STEPS_PARALLEL

    # --------------------------- DEFINE param functions --------------------------------------------
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.itemTsIdReadList = []

    @classmethod
    def worksInStreaming(cls):
        return True

    # --------------------------- INSERT steps functions --------------------------------------------
    def createOutputStep(self, tsId: str):
        with self._lock:
            inSet = self.getInputSet()
            inObj = self.getCurrentItem(inSet, tsId, doLock=False)
            acq = inObj.getAcquisition().clone()
            outputTomos = self._getOutputSetOfTomograms()
            tomo = Tomogram()
            tomo.setTsId(tsId)
            tomo.setFileName(self._getOutTomoFile(tsId))
            tomo.setSamplingRate(inObj.getSamplingRate())
            tomo.setAcquisition(acq)
            tomo.setOrigin()
            if getattr(self, DO_EVEN_ODD, Boolean(False)).get():
                tomo.setHalfMaps([self._getOutTomoFile(tsId, suffix=EVEN), self._getOutTomoFile(tsId, suffix=ODD)])
            outputTomos.append(tomo)
            outputTomos.update(tomo)
            outputTomos.write()
            self._store(outputTomos)
            for outputName in self._possibleOutputs.keys():
                output = getattr(self, outputName, None)
                if output:
                    output.close()

    # --------------------------- INFO functions --------------------------------------------

    # --------------------------- UTILS functions --------------------------------------------
    def rotXTomo(self,
                 tsId: str,
                 suffix: str = '') -> None:
        """Result of the reconstruction must be rotated 90 degrees around the X
        axis to recover the original orientation (due to tomo3d design)"""
        logger.info(cyanStr(f'tsId = {tsId}: rotating the tomogram {suffix.upper()}...'))
        inTomoFile = self._getTmpTomoOutFName(tsId, suffix=suffix)
        outTomoFile = self._getOutTomoFile(tsId, suffix=suffix)

        with mrcfile.mmap(inTomoFile, mode='r', permissive=True) as mrc:
            rotData = np.rot90(mrc.data)

        with mrcfile.mmap(outTomoFile, mode='w+') as mrc:
            mrc.set_data(rotData)

    def readingOutput(self) -> None:
        outTomoSet = getattr(self, self._OUTNAME, None)
        if outTomoSet:
            for item in outTomoSet:
                self.itemTsIdReadList.append(item.getTsId())
            self.info(cyanStr(f'TsIds processed: {self.itemTsIdReadList}'))
        else:
            self.info(cyanStr('No elements have been processed yet'))

    def getTmpFile(self,
                   tsId: str,
                   ext: str = ST_EXT,
                   suf: str = '') -> str:
        baseName = tsId if not suf else f'{tsId}_{suf}'
        ext = ext if ext.startswith('.') else f'.{ext}'
        return self._getTmpPath(f'{baseName}{ext}')

    def getEvenTsTmpFile(self,
                         tsId: str,
                         ext: str = ST_EXT) -> str:
        return self.getTmpFile(tsId, ext, suf=EVEN)

    def getOddTsTmpFile(self,
                        tsId: str,
                        ext: str = ST_EXT) -> str:
        return self.getTmpFile(tsId, ext, suf=ODD)

    def getTsTmpFile(self,
                     tsId: str,
                     ext: str = ST_EXT) -> str:
        return self.getTmpFile(tsId, ext)

    def getWorkingDirName(self, tsId: str) -> str:
        return self._getTmpPath(tsId)

    def _getOutputSetOfTomograms(self) -> SetOfTomograms:
        outTomograms = getattr(self, self._OUTNAME, None)
        if outTomograms:
            outTomograms.enableAppend()
        else:
            inSetPointer = self.getInputSet(pointer=True)
            outTomograms = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            outTomograms.copyInfo(inSetPointer.get())
            outTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{self._OUTNAME: outTomograms})
            self._defineSourceRelation(inSetPointer, outTomograms)
        return outTomograms

    def _getOutTomoFile(self,
                        tsId: str,
                        suffix: str='') -> str:
        return join(self._getTsExtraDir(tsId), f'{tsId}{self._manageSuffix(suffix)}.mrc')

    def _getTmpTomoOutFName(self,
                            tsId: str,
                            suffix: str='') -> str:
        return join(self._getTsTmpDir(tsId), f'{tsId}{self._manageSuffix(suffix)}.mrc')

    @staticmethod
    def _manageSuffix(suffix: str) -> str:
        return f'_{suffix}' if suffix else ''

    def _getTsTmpDir(self, tsId: str) -> str:
        return self._getTmpPath(tsId)

    def _getTsExtraDir(self, tsId: str) -> str:
        return self._getExtraPath(tsId)

    def getCurrentItem(self,
                       inSet: Union[SetOfTomograms, SetOfTiltSeries],
                       tsId: str,
                       doLock: bool = True) -> Union[Tomogram, TiltSeries]:
        if doLock:
            with self._lock:
                return inSet.getItem(Tomogram.TS_ID_FIELD, tsId)
        else:
            return inSet.getItem(Tomogram.TS_ID_FIELD, tsId)

    def getInputSet(self, pointer: bool = False) -> Union[SetOfTomograms, SetOfTiltSeries]:
        # To be defined by the subclasses
        pass
