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
import typing
from enum import Enum
from os.path import join
from typing import Union, List

import mrcfile
import numpy as np
from pyworkflow.object import Set, Boolean, Pointer
from pyworkflow.protocol import STEPS_PARALLEL, IntParam
from pyworkflow.utils import cyanStr, redStr
from tomo.protocols import ProtTomoBase
from pwem.protocols import EMProtocol
from tomo.objects import Tomogram, SetOfTomograms, SetOfTiltSeries, TiltSeries

logger = logging.getLogger(__name__)

# Inputs
IN_TS_SET = 'inputSetOfTiltSeries'
IN_TOMO_SET = 'inputSetTomograms'

# Outputs
OUTPUT_TS_FAILED_NAME = "FailedTiltSeries"
OUTPUT_TOMOS_FAILED_NAME = "FailedTomograms"

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
        self.failedItems = []

    @classmethod
    def worksInStreaming(cls):
        return True

    @staticmethod
    def _insertBinThreadsParam(form):
        form.addParam('binThreads', IntParam,
                      label='threads',
                      default=4,
                      help='Number of threads used by tomo3d each time it is called in the protocol execution. For '
                           'example, if 2 Scipion threads and 3 tomo3d threads are set, the tomograms will be '
                           'processed in groups of 2 at the same time with a call of tomo3d with 3 threads each, so '
                           '6 threads will be used at the same time. Beware the memory of your machine has '
                           'memory enough to load together the number of tomograms specified by Scipion threads.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def createOutputStep(self,
                         inPointer: Pointer,
                         tsId: str):
        try:
            with self._lock:
                inSet = inPointer.get()
                inObj = self.getCurrentItem(inSet, tsId, doLock=False)
                # Set of tomograms
                outputTomos = self.getOutputSetOfTomograms(inPointer)
                # Tomograms
                tomo = Tomogram(tsId=tsId)
                tomo.copyInfo(inObj)
                tomo.setFileName(self._getOutTomoFile(tsId))
                tomo.setOrigin()
                if getattr(self, DO_EVEN_ODD, Boolean(False)).get():
                    tomo.setHalfMaps([self._getOutTomoFile(tsId, suffix=EVEN),
                                      self._getOutTomoFile(tsId, suffix=ODD)])
                outputTomos.append(tomo)
                outputTomos.update(tomo)
                # Data persistence
                outputTomos.write()
                self._store(outputTomos)
                # Close explicitly the outputs (for streaming)
                for outputName in self._possibleOutputs.keys():
                    output = getattr(self, outputName, None)
                    if output:
                        output.close()
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with '
                                f'exception {e}. Skipping... '))

    def closeOutputSetsStep(self, attrib: Union[List[str], str]):
        self._closeOutputSet()
        attribList = [attrib] if type(attrib) is str else attrib
        failedOutputList = []
        for attr in attribList:
            outTsSet = getattr(self, attr, None)
            if not outTsSet or (outTsSet and len(outTsSet) == 0):
                failedOutputList.append(attr)
        if failedOutputList:
            raise Exception(f'No output/s {failedOutputList} were generated. Please check the '
                            f'Output Log > run.stdout and run.stderr')

    # --------------------------- OUTPUTS functions --------------------------------------------
    def getOutputSetOfTomograms(self,
                                inputPtr: Pointer) -> SetOfTomograms:
        inputSet = inputPtr.get()
        attribName = outputTomo3dObjects.tomograms.name
        outputSet = getattr(self, attribName, None)

        if outputSet:
            outputSet.enableAppend()
        else:
            outputSet = self._createSetOfTomograms()
            outputSet.copyInfo(inputSet)
            outputSet.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{attribName: outputSet})
            self._defineSourceRelation(inputPtr, outputSet)

        return outputSet

    def getOutputFailedSet(self,
                           inputPtr: Pointer,
                           inputsAreTs: bool = True) -> Union[SetOfTiltSeries, SetOfTomograms]:
        """ Create output set for failed TS or tomograms. """
        inputSet = inputPtr.get()
        if inputsAreTs:
            failedTs = getattr(self, OUTPUT_TS_FAILED_NAME, None)

            if failedTs:
                failedTs.enableAppend()
            else:
                logger.info(cyanStr('Create the set of failed TS'))
                failedTs = self._createSetOfTiltSeries(suffix='Failed')
                failedTs.copyInfo(inputSet)
                failedTs.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**{OUTPUT_TS_FAILED_NAME: failedTs})
                self._defineSourceRelation(inputPtr, failedTs)

            return failedTs

        else:
            failedTomos = getattr(self, OUTPUT_TOMOS_FAILED_NAME, None)
            if failedTomos:
                failedTomos.enableAppend()
            else:
                logger.info(cyanStr('Create the set of failed tomograms'))
                failedTomos = self._createSetOfTomograms(suffix='Failed')
                failedTomos.copyInfo(inputSet)
                failedTomos.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**{OUTPUT_TOMOS_FAILED_NAME: failedTomos})
                self._defineSourceRelation(inputPtr, failedTomos)

            return failedTomos

    def addToOutFailedSet(self,
                          tsId: str,
                          inputsAreTs: bool = True) -> None:
        """ Just copy input item to the failed output set. """
        logger.info(cyanStr(f'Failed TS ---> {tsId}'))
        try:
            inputSet = self.getInputTsSet(pointer=True) if inputsAreTs else self.getInputTomoSet(pointer=True)
            output = self.getOutputFailedSet(inputSet, inputsAreTs=inputsAreTs)
            with self._lock:
                item = self.getCurrentTs(tsId) if inputsAreTs else self.getCurrentTomo(tsId)
                newItem = item.clone()
                newItem.copyInfo(item)
                output.append(newItem)

                if isinstance(item, TiltSeries):
                    newItem.copyItems(item)
                    newItem.write()

                output.update(newItem)
                output.write()
                self._store(output)
                # Close explicitly the outputs (for streaming)
                output.close()
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the failed output with '
                                f'exception {e}. Skipping... '))

    # --------------------------- INFO functions --------------------------------------------

    # --------------------------- UTILS functions --------------------------------------------
    def getInputTsSet(self, pointer: bool = False) -> typing.Union[Pointer, SetOfTiltSeries]:
        tsSetPointer = getattr(self, IN_TS_SET)
        return tsSetPointer if pointer else tsSetPointer.get()

    def getCurrentTs(self, tsId: str) -> TiltSeries:
        return self.getInputTsSet().getItem(TiltSeries.TS_ID_FIELD, tsId)

    def getInputTomoSet(self, pointer: bool = False) -> Union[Pointer, SetOfTomograms]:
        tomoSetPointer = getattr(self, IN_TOMO_SET)
        return tomoSetPointer if pointer else tomoSetPointer.get()

    def getCurrentTomo(self, tsId: str) -> Tomogram:
        return self.getInputTomoSet().getItem(Tomogram.TS_ID_FIELD, tsId)

    def getCurrentItem(self,
                       inSet: Union[SetOfTomograms, SetOfTiltSeries],
                       tsId: str,
                       doLock: bool = True) -> Union[Tomogram, TiltSeries]:
        if doLock:
            with self._lock:
                return inSet.getItem(Tomogram.TS_ID_FIELD, tsId)
        else:
            return inSet.getItem(Tomogram.TS_ID_FIELD, tsId)

    def readingOutput(self) -> None:
        outTomoSet = getattr(self, self._OUTNAME, None)
        if outTomoSet:
            for item in outTomoSet:
                self.itemTsIdReadList.append(item.getTsId())
            self.info(cyanStr(f'TsIds processed: {self.itemTsIdReadList}'))
        else:
            self.info(cyanStr('No elements have been processed yet'))

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

