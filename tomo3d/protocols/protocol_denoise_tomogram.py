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
import time
import typing
from pyworkflow.object import Pointer
from pyworkflow.protocol import ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, makePath
from tomo.objects import SetOfTomograms
from tomo3d import Plugin
from pyworkflow.protocol.params import IntParam, EnumParam, LEVEL_ADVANCED, FloatParam, PointerParam, GT
from tomo3d.protocols.protocol_base import ProtBaseTomo3d

logger = logging.getLogger(__name__)
IN_TOMOS = 'inputSetTomograms'

#Denoising methods
DENOISE_EED = 0
DENOISE_BF = 1


class ProtTomo3dProtDenoiseTomogram(ProtBaseTomo3d, ProtStreamingBase):
    """ Denoises sets of tomograms using methods described in https://sites.google.com/site/3demimageprocessing/
    Two methods are available: \n
    _TomoAND (also known as TomoEED)_ is an optimized program for denoising tomographic volumes with
    the Anisotropic Nonlinear Diffusion (AND) method, using the EED (Edge-Enhancing Diffusion) mode.
    This mode manages to reduce noise in the volume with good abilites to preserve and enhance the structures.\n
    _TomoBflow_: The program TOMOBFLOW (pronounced as tomo-be-flow) is intended for noise filtering with preservation
     of biologically relevant information. It is an efficient implementation of the Beltrami flow, a nonlinear
      filtering method that locally tunes the strength of the smoothing according to an edge indicator based
      on geometry properties. TOMOBFLOW is equipped with the power of diffusion-based filtering methods, with
      the important advantage that it does not require complicated parameter tuning for successful denoising
      of datasets. Furthermore, the program has been optimized to reduce the computational demands, specially
      in terms of memory requirements.
    """
    _label = 'denoise tomogram'

    def __init__(self, **args):
        super().__init__(**args)

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMOS, PointerParam, pointerClass='SetOfTomograms',
                      label='Set Of Tomograms',
                      help='Set of tomograms that will be denoised.')
        form.addParam('method', EnumParam,
                      choices=['Edge Enhancing Diffusion (EED)', 'BFlow'],
                      default=DENOISE_EED,
                      label='Denoising method',
                      help='Denoising method to use')
        form.addParam('SigmaGaussian', FloatParam, default=0.5,
                      label='Sigma Gaussian Filter',
                      help='This option allows the user to enter the standard deviation '
                           'for an initial Gaussian filtering of the input tomogram. '
                           'This initial filtering is intended to remove the fine-grain noise '
                           'in order to regularize the computation of the gradients and '
                           'make them less sensitive to noise. A value in the range [0.5,1.0] '
                           'is advisable. Higher values may blur the structures, which '
                           'would thus spoil the effects of the subsequent Beltrami flow. '
                           'By default a value of 0.5 is assumed. If no initial Gaussian '
                           'filtering is wanted, it must be switched off with "-g 0".')

        form.addParam('nIter', IntParam, default=10,
                      label='Number of Iterations',
                      condition='method == %i' % DENOISE_EED,
                      help='Number of diffusion iterations.\n\n'
                           '_Edge Enhancing Diffusion (EED)_: Number of diffusion iterations.'
                           ' A number of iterations around [10,60] are OK, though it '
                           ' strongly depends on the parameter K (lambda). If K is set '
                           ' to time-varying average gradient in the whole tomogram '
                           ' (i.e. the default behaviour), the number of iterations should '
                           ' be relatively low (10–20). The more iterations the stronger the '
                           ' smoothing. If too many iterations are used, some interesting '
                           ' structural features may turn out to be blurred. By default, '
                           ' 10 iterations are assumed. This parameter is strongly related'
                           ' to the time step \n')

        form.addParam('nIterBflow', IntParam,
                      default=70,
                      validators=[GT(0)],
                      label='Number of Iterations',
                      condition='method == %i' % DENOISE_BF,
                      help='Number of diffusion iterations.\n\n'
                           '_BFlow Iterations_: A number of iterations around [50,250] yields '
                           ' fairly good results. The more iterations the stronger the smoothing. '
                           ' If too many iterations are used, some interesting structural '
                           ' features may turn out to be blurred. By default, 70 iterations '
                           ' are assumed. This parameter is strongly related to the time '
                           ' step.\n')
        #
        form.addParam('Lambda', FloatParam, default=-1.0,
                      condition='method == %i' % DENOISE_EED,
                      label='Lambda',
                      help='This option allows the user to specify a value for the '
                           ' parameter K (lambda) of Anisotropic Nonlinear Diffusion. '
                           ' This value will be used throughout the denoising process. '
                           ' That is, it remains constant for all iterations. By default, '
                           ' the program uses a time-varying K set as the average '
                           ' gradient computed from the whole tomogram.',
                      expertLevel=LEVEL_ADVANCED)

        form.addParam('TimeStep', FloatParam, default=0.1,
                      label='Time Step',
                      condition='method == %i' % DENOISE_EED,
                      help='This option allows the user to change the time step. The larger '
                           ' the time step, the lower the number of iterations needed. '
                           ' The default value is 0.1, which is the standard value. However, '
                           ' the maximum value for the sake of numerical stability that still '
                           ' works for explicit discretization of partial differential '
                           ' equations is 0.15. So, the user might want to set this parameter '
                           ' to 0.15, and decrease the number of iterations accordingly.',
                      expertLevel=LEVEL_ADVANCED)

        form.addParam('TimeStepBflow', FloatParam, default=0.15,
                      label='Time Step',
                      condition='method == %i' % DENOISE_BF,
                      help='This option allows the user to change the time step. The default '
                           ' value is 0.15, which is the maximum value for the sake of '
                           ' numerical stability. Recommended values should be in the '
                           ' range [0.1,0.15]. The larger the time step, the lower the '
                           ' number of iterations needed.',
                      expertLevel=LEVEL_ADVANCED)
        self._insertBinThreadsParam(form)
        form.addParallelSection(threads=2, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def stepsGeneratorStep(self) -> None:
        closeSetStepDeps = []
        inTomoSet = self.getInputSet()
        self.readingOutput()

        while True:
            listInTsIds = inTomoSet.getTSIds()
            if not inTomoSet.isStreamOpen() and self.itemTsIdReadList == listInTsIds:
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self._closeOutputSet,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break
            for tomo in inTomoSet.iterItems():
                tsId = tomo.getTsId()
                if tsId not in self.itemTsIdReadList:
                    cInputId = self._insertFunctionStep(self.denoiseTomogramStep, tsId,
                                                        prerequisites=[],
                                                        needsGPU=False)
                    cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                                      prerequisites=cInputId,
                                                      needsGPU=False)
                    closeSetStepDeps.append(cOutId)
                    logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                    self.itemTsIdReadList.append(tsId)
            time.sleep(10)
            if inTomoSet.isStreamOpen():
                with self._lock:
                    inTomoSet.loadAllProperties()  # refresh status for the streaming

    # --------------------------- STEPS functions --------------------------------------------
    def denoiseTomogramStep(self, tsId: str):
        makePath(self._getTsExtraDir(tsId))
        if self.method.get() == DENOISE_EED:
            logger.info(cyanStr(f'tsId = {tsId}: denoising by Edge Enhancing Diffusion...'))
            self.runEED(tsId)
        else:  # self.method.get() == DENOISE_BF:
            logger.info(cyanStr(f'tsId = {tsId}: denoising by BFlow...'))
            self.runBflow(tsId)

    # --------------------------- INFO functions --------------------------------------------
    def _citations(self):
        return ['Fernandez2018_tomoeed', 'Fernandez2009_tomobflow']

    # --------------------------- UTILS functions --------------------------------------------
    def getInputSet(self, pointer: bool = False) -> typing.Union[Pointer, SetOfTomograms]:
        tomoSetPointer = getattr(self, IN_TOMOS)
        return tomoSetPointer if pointer else tomoSetPointer.get()

    def runBflow(self, tsId: str):
        """Denoises de tomogram using the BFlow method"""
        params = self.getCommonParamsCmd(tsId)
        params.append(f'-i {self.nIterBflow.get()}')
        args = ' '.join(params)
        self.runJob(Plugin.getTomoBFlowProgram(), args)

    def runEED(self, tsId: str):
        """Denoises de tomogram using the EED method"""
        params = self.getCommonParamsCmd(tsId)
        params.append(f'-i {self.nIter.get()}')
        if self.Lambda.get() >= 0:
            params.append(f'-k {self.Lambda.get()}')
        args = ' '.join(params)
        self.runJob(Plugin.getTomoEEDProgram(), args)

    def getCommonParamsCmd(self, tsId: str) -> typing.List[str]:
        tomo = self.getCurrentItem(self.getInputSet(), tsId)
        inTomoFile = tomo.getFileName()
        outTomoFile = self._getOutTomoFile(tsId)
        params = [
            f'{inTomoFile}',
            f'{outTomoFile}'
            f' -g {self.SigmaGaussian.get()}',
            f'-s {self.TimeStep.get()}',
            f'-t {self.binThreads.get()}'
        ]
        return params


