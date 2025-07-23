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
import logging
import time
import typing

from pyworkflow.object import Pointer, Set
from pyworkflow.protocol import ProtStreamingBase
from tomo.objects import SetOfTiltSeries, TiltSeries
from tomo3d.protocols.protocol_base import ProtBaseTomo3d, EVEN, ODD, DO_EVEN_ODD, RAWTLT_EXT
from pyworkflow.utils import makePath, Message, cyanStr, redStr
from tomo3d import Plugin
from pyworkflow.protocol.params import IntParam, EnumParam, FloatParam, LEVEL_ADVANCED, BooleanParam, PointerParam, GE, \
    GT

logger = logging.getLogger(__name__)
IN_TS_SET = 'inputSetOfTiltSeries'
OUTPUT_TS_FAILED_NAME = "FailedTiltSeries"

# Reconstruction methods
WBP = 0
SIRT = 1


class ProtTomo3dReconstrucTomo(ProtBaseTomo3d, ProtStreamingBase):
    """ Reconstruct tomograms using TOMO3D from
    Software from: https://sites.google.com/site/3demimageprocessing/

    Tomo3D implements a multithreaded vectorized approach to tomographic reconstruction
    that takes full advantage of the power in modern multicore computers. Full resolution
    tomograms are generated at high speed on standard computers with no special system
    requirements. Tomo3D has the most common reconstruction methods implemented,
    namely WBP and SIRT. It proves to be competitive with current GPU solutions in terms
    of processing time, in the order of a few seconds with WBP or minutes with SIRT.
    """
    _label = 'reconstruct tomogram'

    def __init__(self, **args):
        super().__init__(**args)

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineInputParams(form)
        form.addParam('method', EnumParam,
                      choices=['WBP (Fast)', 'SIRT (Slow)'], default=WBP,
                      label='Reconstruction method',
                      help='Reconstrution method to use:\n'
                           '_WBP_: Also called filtered back projection provides an exact'
                           'solution for the reconstruction problem. Unfortunately, the '
                           'contrast of the reconstructed tomogram is poor. This solution should'
                           'be chosen for subtomogram averaging approaches\n'
                           '_SIRT_: It is an algebraic iterative algorithm that aims to maximize'
                           'the compabitity between the images and the reconstructed tomogram. This'
                           'solution is more suitable for celullar environment, but not for '
                           'subtomogram averaging')
        form.addParam('nIterations', IntParam,
                      default=30,
                      validators=[GT(0)],
                      condition='method==%i' % SIRT,
                      label='Number of Iterations (SIRT)',
                      help='It sets the number of iterations for SIRT. By default,'
                           ' a number of 30 iterations is used.')
        self._defineSetShapeParams(form)
        form.addParam('Hamming', FloatParam, default=0.0,
                      condition='method==%i' % WBP,
                      expertLevel=LEVEL_ADVANCED,
                      label='Hamming filter frequency',
                      help=' By default, in WBP the program applies a Ramp filter together '
                           ' with a Hamming filter. The latter allows attenuation of the '
                           ' high frequency components, which are typically very noisy. '
                           ' This parameter allows setting up the starting frequency from '
                           ' which the Hamming filter is to be applied. The values '
                           ' must range in [0,0.5], where 0.5 represents the Nyquist frequency.'
                           ' By default, the program uses a starting frequency of 0, i.e. a '
                           ' Hamming is applied over the whole frequency range. A value of 0.5 '
                           ' would turn off the Hamming filter. Values in-between would '
                           ' preserve low-frequency components and modulate the contribution '
                           ' of the other frequency components')
        self._insertBinThreadsParam(form)
        form.addParallelSection(threads=2, mpi=0)

    @staticmethod
    def _defineInputParams(form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TS_SET, PointerParam,
                      important=True,
                      pointerClass='SetOfTiltSeries',
                      label='Tilt series',
                      help='Tilt Series to reconstruct the tomograms. Ideally these tilt series'
                            'should contain alignment information.')
        form.addParam(DO_EVEN_ODD, BooleanParam,
                      label='Reconstruct the even/odd tomograms?',
                      default=False)

    @staticmethod
    def _defineSetShapeParams(form):
        form.addParam('setShape', BooleanParam,
                      default=True,
                      label='Set manual tomogram shape',
                      display=BooleanParam,
                      help='By default the shape of the tomogram is defined by the tilt series shape,'
                           'but the user can manually set it here')

        group = form.addGroup('Tomogram dimensions', condition='setShape')

        group.addParam('height', IntParam,
                       default=300,
                       validators=[GE(0)],
                       label='Thickness (px)',
                       help='By default, the height (i.e. thickness) of the reconstructed tomogram '
                            ' is equal to the X dimension of the images in the tilt-series. '
                            ' This option allows the user to properly set the height of the '
                            ' tomogram to a value that better fits the thickness of the '
                            ' specimen under study. This also allows reduction of the processing '
                            ' time. The thickness value has to be an even number. This is '
                            ' necessary to exploit projection symmetry (see Agulleiro et al. '
                            ' J.Struct.Biol. 170:570–575, 2010).')

        group.addParam('width', IntParam,
                       expertLevel=LEVEL_ADVANCED,
                       default=0,
                       validators=[GE(0)],
                       label='Width (px)',
                       help='By default, the width of the reconstructed tomogram is equal to '
                            ' the X dimension of the images in the tilt-series. This option '
                            ' allows the user to set the width of the tomogram to a smaller '
                            ' value, which helps to better focus the reconstruction to an area '
                            ' of interest and further reduce the processing time. Regardless '
                            ' of the width value, the tomogram is always reconstructed around '
                            ' the tilt axis (i.e. the centre of the images in the tilt-series).')

        line = group.addLine('Slices in the Y-axis',
                             expertLevel=LEVEL_ADVANCED,
                             help=' By default, the output tomogram has as many slices along '
                                  ' the tilt axis as the Y dimension of the images in the tilt-series. '
                                  ' This option allows the user to specify the set of slices to be '
                                  ' reconstructed. Thus, this option helps to better focus the '
                                  ' reconstruction to an area of interest and further reduce the '
                                  ' processing time. The user has to specify two values y1,y2 numbered '
                                  ' from 0 (i.e. in the range [0,Ny-1], with Ny being the Y dimension '
                                  ' of the images), which represent the indices of the first and last '
                                  ' slices to be included in the output tomogram.')

        line.addParam('iniSlice', IntParam, default=0, validators=[GE(0)], label='Initial')
        line.addParam('finSlice', IntParam, default=0, validators=[GE(0)], label='Final')

    # --------------------------- INSERT steps functions --------------------------------------------
    def stepsGeneratorStep(self) -> None:
        closeSetStepDeps = []
        inTsSet = self.getInputSet()
        self.readingOutput()

        while True:
            listInTsIds = inTsSet.getTSIds()
            if not inTsSet.isStreamOpen() and self.itemTsIdReadList == listInTsIds:
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self._closeOutputSet,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break
            for ts in inTsSet.iterItems():
                tsId = ts.getTsId()
                if tsId not in self.itemTsIdReadList and ts.getSize() > 0:  # Avoid processing empty TS (before the Tis are added)
                    cInputId = self._insertFunctionStep(self.convertInputStep, tsId,
                                                        prerequisites=[],
                                                        needsGPU=False)
                    recId = self._insertFunctionStep(self.reconstructTomogramStep, tsId,
                                                     prerequisites=cInputId,
                                                     needsGPU=False)
                    cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                                      prerequisites=recId,
                                                      needsGPU=False)
                    closeSetStepDeps.append(cOutId)
                    logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                    self.itemTsIdReadList.append(tsId)
            time.sleep(10)
            if inTsSet.isStreamOpen():
                with self._lock:
                    inTsSet.loadAllProperties()  # refresh status for the streaming

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self, tsId: str):
        logger.info(cyanStr(f'tsId = {tsId}: converting the inputs...'))
        makePath(self.getWorkingDirName(tsId), self._getTsExtraDir(tsId))
        tsTmpFile = self.getTsTmpFile(tsId)
        tltTmpFile = self.getTsTmpFile(tsId, ext=RAWTLT_EXT)
        ts = self.getCurrentItem(self.getInputSet(), tsId)
        rotationAngle = ts.getAcquisition().getTiltAxisAngle()
        # Check if rotation angle is greater than 45º. If so,
        # swap x and y dimensions to adapt output image sizes to
        # the final sample disposition.
        swapXY = True if 45 < abs(rotationAngle) < 135 else False
        if self.doEvenOdd.get():
            tsEvenTmpFile = self.getEvenTsTmpFile(tsId)
            tsOddTmpFile = self.getOddTsTmpFile(tsId)
            ts.applyTransformToAll(tsTmpFile,
                                   swapXY=swapXY,
                                   outFileNamesEvenOdd=[tsEvenTmpFile, tsOddTmpFile])
        else:
            ts.applyTransform(tsTmpFile,
                              swapXY=swapXY)
        ts.generateTltFile(tltTmpFile)

    def reconstructTomogramStep(self, tsId: str):
        try:
            params = ''
            if self.method.get() == SIRT:
                params += f' -S -l {self.nIterations.get()}'
            if self.setShape.get():
                if self.width.get() != 0:
                    params += f' -x {self.width.get()}'
                if self.finSlice.get() != 0:
                    params += f' -y {self.iniSlice.get()},{self.finSlice.get()}'
                if self.height.get() != 0:
                    params += f' -z {self.height.get()}'

            tomoRecInfoDict = {'': self.getTsTmpFile(tsId)}  # {suffix: fileName}
            tltTmpFile = self.getTsTmpFile(tsId, ext=RAWTLT_EXT)
            if self.doEvenOdd.get():
                tomoRecInfoDict[EVEN] = self.getEvenTsTmpFile(tsId)
                tomoRecInfoDict[ODD] = self.getOddTsTmpFile(tsId)

            for suffix, inFile in tomoRecInfoDict.items():
                logger.info(cyanStr(f'tsId = {tsId}: reconstructing the tomogram {suffix.upper()}...'))
                outFile = self._getTmpTomoOutFName(tsId, suffix=suffix)
                args = f'-i {inFile} -a {tltTmpFile} -o {outFile} -t {self.binThreads.get()}'
                args += params
                self.runJob(Plugin.getTomo3dProgram(), args)
                self.rotXTomo(tsId, suffix=suffix)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'Tomo3d tomogram reconstruction execution failed for tsId {tsId} -> {e}'))

    def createOutputStep(self, tsId: str):
        if tsId in self.failedItems:
            with self._lock:
                inTs = self.getInputSet().getItem(TiltSeries.TS_ID_FIELD, tsId)
                self.createOutputFailedSet(inTs)
                failedTs = getattr(self, OUTPUT_TS_FAILED_NAME, None)
                if failedTs:
                    failedTs.close()
        else:
            super().createOutputStep(tsId)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errorMsg = []
        if self.height.get() % 2 == 1:
            errorMsg.append('The thickness must be an even number')
        if self.doEvenOdd.get() and not self.inputSetOfTiltSeries.get().hasOddEven():
            errorMsg.append('The even/odd tomograms cannot be reconstructed as no even/odd tilt-series are found '
                            'in the metadata of the introduced tilt-series.')

        return errorMsg

    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Fernandez2010_tomo3d', 'Fernandez2015_tomo3d']

# --------------------------- UTILS functions --------------------------------------------
    def getInputSet(self, pointer: bool = False) -> typing.Union[Pointer, SetOfTiltSeries]:
        tsSetPointer = getattr(self, IN_TS_SET)
        return tsSetPointer if pointer else tsSetPointer.get()

    def createOutputFailedSet(self, item):
        """ Just copy input item to the failed output set. """
        logger.info(f'Failed TS ---> {item.getTsId()}')
        inputSetPointer = self.getInputSet(pointer=True)
        output = self.getOutputFailedSet(inputSetPointer)
        newItem = item.clone()
        newItem.copyInfo(item)
        output.append(newItem)

        if isinstance(item, TiltSeries):
            newItem.copyItems(item)
            newItem.write(properties=False)

        output.update(newItem)
        output.write()
        self._store(output)

    def getOutputFailedSet(self, inputPtr: Pointer):
        """ Create output set for failed TS or tomograms. """
        inputSet = inputPtr.get()
        if isinstance(inputSet, SetOfTiltSeries):
            failedTs = getattr(self, OUTPUT_TS_FAILED_NAME, None)

            if failedTs:
                failedTs.enableAppend()
            else:
                logger.info('Create the set of failed TS')
                failedTs = SetOfTiltSeries.create(self._getPath(), template='tiltseries', suffix='Failed')
                failedTs.copyInfo(inputSet)
                failedTs.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**{OUTPUT_TS_FAILED_NAME: failedTs})
                self._defineSourceRelation(inputPtr, failedTs)

            return failedTs

