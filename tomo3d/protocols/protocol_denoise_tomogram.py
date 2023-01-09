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
from enum import Enum
from os.path import basename

from pyworkflow import BETA
from tomo.protocols import ProtTomoBase
from tomo3d import Plugin

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import IntParam, EnumParam, LEVEL_ADVANCED, FloatParam, PointerParam

from tomo.objects import Tomogram, SetOfTomograms

DENOISE_EED = 0
DENOISE_BF = 1


class outputDenoiseObjects(Enum):
    tomograms = SetOfTomograms


class ProtJjsoftProtDenoiseTomogram(EMProtocol, ProtTomoBase):
    """ Denoises sets of tomograms using methods described in https://sites.google.com/site/3demimageprocessing/
    Returns the set of denoised tomograms
    """
    _label = 'denoise tomogram'
    _devStatus = BETA
    _possibleOutputs = outputDenoiseObjects

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label='Input')
        form.addParam('inputSetTomograms', PointerParam, pointerClass='SetOfTomograms',
                      label='Set Of Tomograms',
                      help='Select one set of tomograms')
        form.addParam('method', EnumParam,
                      choices=['Edge Enhancing Diffusion (EED)', 'BFlow'],
                      default=DENOISE_EED,
                      label='Denoising method',
                      help='Denoising method to use')
        form.addSection(label='Parameters')
        form.addParam('SigmaGaussian', FloatParam, default=0.5,
                      label='Sigma Gaussian Filter',
                      help='Sigma for initial gaussian filtering.')
        form.addParam('nIter', IntParam, default=10,
                      label='Number of Iterations',
                      help='Number of Iterations of denoising.')
        #
        form.addParam('Lambda', FloatParam, default=-1.0,
                      condition = 'method==0',
                      label='Lambda (EED)',
                      help='Lambda threshold for gaussian filtering (Default (Negative): time-varying value estimated)',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('TimeStep', FloatParam, default=0.1,
                      label='Time Step',
                      help='Time Step for Iterations (max 0.15)',
                      expertLevel=LEVEL_ADVANCED)
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert every step of the protocol"""
        inputTomos = self.inputSetTomograms.get()
        self.outputFiles = []
        pre = []
        for tomo in inputTomos.iterItems():
            stepId = self._insertFunctionStep(self.denoiseTomogramStep, tomo.getFileName())
            pre.append(stepId)

        self._insertFunctionStep(self.createOutputStep, prerequisites=pre)

    # --------------------------- STEPS functions --------------------------------------------
    def denoiseTomogramStep(self, inp_tomo_path):
        # We start preparing writing those elements we're using as input to keep them untouched
        if self.method.get() == DENOISE_EED:
            print('Denoising by Edge Enhancing Diffusion')
            # call EED
            out_tomo_path = self.call_EED(inp_tomo_path)

        else:  # self.method.get() == DENOISE_BF:
            print('Denoising by BFlow')
            # call BFlow
            out_tomo_path = self.call_BFlow(inp_tomo_path)

        self.outputFiles.append(out_tomo_path)

    def createOutputStep(self):
        inputTomos = self.inputSetTomograms.get()
        outputTomos = self._createSetOfTomograms()
        outputTomos.copyInfo(inputTomos)

        for i, inp_tomo in enumerate(inputTomos):
            tomo_path = self.outputFiles[i]
            tomo = Tomogram()
            tomo.copyInfo(inp_tomo)
            tomo.setLocation(tomo_path)
            outputTomos.append(tomo)

        self._defineOutputs(**{outputDenoiseObjects.tomograms.name: outputTomos})
        self._defineSourceRelation(self.inputSetTomograms, outputTomos)

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        pass

    def _methods(self):
        pass

    def _citations(self):
        return ['Fernandez2018_tomoeed', 'Fernandez2009_tomobflow']

    # --------------------------- UTILS functions --------------------------------------------
    def call_BFlow(self, inp_tomo_path):
        """Denoises de tomogram using the AND method"""
        params = '-g {} -i {} -s {} -t {}'.format(self.SigmaGaussian.get(), self.nIter.get(),
                                                  self.TimeStep.get(), self.numberOfThreads)
        out_tomo_path = self._getExtraPath(basename(inp_tomo_path))
        args = '{} {} {}'.format(params, inp_tomo_path, out_tomo_path)
        self.runJob(Plugin.getTomoBFlowProgram(), args)
        return out_tomo_path

    def call_EED(self, inp_tomo_path):
        """Denoises de tomogram using the AND method"""
        if self.Lambda.get() < 0:
            params = '-g {} -i {} -s {}'.format(self.SigmaGaussian.get(), self.nIter.get(), self.TimeStep.get())
        else:
            params = '-g {} -i {} -s {} -k {}'.format(self.SigmaGaussian.get(), self.nIter.get(),
                                                      self.TimeStep.get(),self.Lambda.get())
        out_tomo_path = self._getExtraPath(basename(inp_tomo_path))
        args = '{} {} {}'.format(params, inp_tomo_path, out_tomo_path)
        self.runJob(Plugin.getTomoEEDProgram(), args)
        return out_tomo_path


