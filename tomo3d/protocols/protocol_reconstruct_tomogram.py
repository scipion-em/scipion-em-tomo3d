# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
# *****************************************************************************
import logging
import time
from collections import Counter

from pyworkflow.protocol import ProtStreamingBase
from tomo.objects import TiltSeries
from tomo3d.protocols.protocol_base import ProtBaseTomo3d, EVEN, ODD, DO_EVEN_ODD, RAWTLT_EXT, outputTomo3dObjects, \
    IN_TS_SET
from pyworkflow.utils import makePath, Message, cyanStr, redStr
from tomo3d import Plugin
from pyworkflow.protocol.params import IntParam, EnumParam, FloatParam, LEVEL_ADVANCED, BooleanParam, PointerParam, GE, \
    GT

logger = logging.getLogger(__name__)

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
    program = Plugin.getTomo3dProgram()

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
        inTsSet = self.getInputTsSet()
        self.readingOutput()

        while True:
            with self._lock:
                inTsIds = set(inTsSet.getTSIds())

            # In the if statement below, Counter is used because in the tsId comparison the order doesn’t matter
            # but duplicates do. With a direct comparison, the closing step may not be inserted because of the order:
            # ['ts_a', 'ts_b'] != ['ts_b', 'ts_a'], but they are the same with Counter.
            if not inTsSet.isStreamOpen() and Counter(self.itemTsIdReadList) == Counter(inTsIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         self._OUTNAME,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            nonProcessedTsIds = inTsIds - set(self.itemTsIdReadList)
            tsToProcessDict = {tsId: ts.clone() for ts in inTsSet.iterItems()
                               if (tsId := ts.getTsId()) in nonProcessedTsIds  # Only not processed tsIds
                               and ts.getSize() > 0}  # Avoid processing empty TS
            for tsId, ts in tsToProcessDict.items():
                cInputId = self._insertFunctionStep(self.convertInputStep, ts,
                                                    prerequisites=[],
                                                    needsGPU=False)
                recId = self._insertFunctionStep(self.reconstructTomogramStep, tsId,
                                                 prerequisites=cInputId,
                                                 needsGPU=False)
                cOutId = self._insertFunctionStep(self.createOutStep, ts,
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
    def convertInputStep(self, ts: TiltSeries):
        tsId = ts.getTsId()
        try:
            logger.info(cyanStr(f'tsId = {tsId}: converting the inputs...'))
            makePath(self.getWorkingDirName(tsId), self._getTsExtraDir(tsId))
            tsTmpFile = self.getTsTmpFile(tsId)
            tltTmpFile = self.getTsTmpFile(tsId, ext=RAWTLT_EXT)
            if self.doEvenOdd.get():
                tsEvenTmpFile = self.getEvenTsTmpFile(tsId)
                tsOddTmpFile = self.getOddTsTmpFile(tsId)
                ts.applyTransformToAll(tsTmpFile,
                                       outFileNamesEvenOdd=[tsEvenTmpFile, tsOddTmpFile])
            else:
                ts.applyTransform(tsTmpFile)
            ts.generateTltFile(tltTmpFile)
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed with the exception -> {e}'))

    def reconstructTomogramStep(self, tsId: str):
        if tsId not in self.failedItems:
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
                    self.runJob(self.program, args)
                    self.rotXTomo(tsId, suffix=suffix)
            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {self.program} execution failed'
                                    f' with the exception -> {e}'))

    def createOutStep(self, ts: TiltSeries):
        tsId = ts.getTsId()
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
        else:
            inTsPointer = self.getInputTsSet(pointer=True)
            super().createOutputStep(inTsPointer, ts)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errorMsg = []
        if self.height.get() % 2 == 1:
            errorMsg.append('The thickness must be an even number')
        if self.doEvenOdd.get() and not self.getInputTsSet().hasOddEven():
            errorMsg.append('The even/odd tomograms cannot be reconstructed as no even/odd tilt-series are found '
                            'in the metadata of the introduced tilt-series.')

        return errorMsg

    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Fernandez2010_tomo3d', 'Fernandez2015_tomo3d']

# --------------------------- UTILS functions --------------------------------------------



