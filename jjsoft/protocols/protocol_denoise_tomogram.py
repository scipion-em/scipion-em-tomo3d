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

from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import IntParam, EnumParam, LEVEL_ADVANCED, FloatParam, BooleanParam, PointerParam

from tomo.objects import Tomogram, SetOfTomograms

class JjsoftProtDenoiseTomogram(EMProtocol):
    """ Remove particles noise by filtering them.
    This filtering process is based on a projection over a basis created
    from some averages (extracted from classes). This filtering is not
    intended for processing particles. The huge filtering they will be
    passed through is known to remove part of the signal with the noise.
    However this is a good method for clearly see which particle are we
    going to process before it's done.
    """
    _label = 'denoise tomogram'

    DENOISE_EED = 0
    DENOISE_BF = 1

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
                      choices=['Edge Enhancing Diffusion (EED)','BFlow'], default=0,
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
        iter = 1
        for tomo in inputTomos.iterItems():
            stepId= self._insertFunctionStep('denoiseTomogramStep', tomo.getFileName(), iter)
            pre.append(stepId)
            iter += 1

        self._insertFunctionStep('createOutputStep', prerequisites=pre)

    # --------------------------- STEPS functions --------------------------------------------
    def denoiseTomogramStep(self, inp_tomo_path, iter):
        # We start preparing writing those elements we're using as input to keep them untouched
        if self.method.get() == 0:
            print('Denoising by Edge Enhancing Diffusion')
            #Call EED
            out_tomo_path = self.call_EED(inp_tomo_path)

        elif  self.method.get() == 1:
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
            tomo.setLocation(tomo_path)
            tomo.setOrigin(inp_tomo.getOrigin())
            tomo.setAcquisition(inp_tomo.getAcquisition())
            outputTomos.append(tomo)

        self._defineOutputs(outputTomograms=outputTomos)
        self.outputTomograms=outputTomos
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
        return ['Fernandez2018','Fernandez2009']

    # --------------------------- UTILS functions --------------------------------------------
    def call_BFlow(self, inp_tomo_path):
        '''Denoises de tomogram using the AND method'''
        params = '-g {} -i {} -s {} -t {}'.format(self.SigmaGaussian.get(), self.nIter.get(),
                                                  self.TimeStep.get(), self.numberOfThreads)
        out_tomo_path = self._getExtraPath('denoisedBflow_'+inp_tomo_path.split('/')[-1])
        args = '{} {} {}'.format(params, inp_tomo_path, out_tomo_path)
        self.runJob('tomobflow', args)
        return out_tomo_path

    def call_EED(self, inp_tomo_path):
        '''Denoises de tomogram using the AND method'''
        if self.Lambda.get() < 0:
            params = '-g {} -i {} -s {}'.format(self.SigmaGaussian.get(), self.nIter.get(), self.TimeStep.get())
        else:
            params = '-g {} -i {} -s {} -k {}'.format(self.SigmaGaussian.get(), self.nIter.get(),
                                                      self.TimeStep.get(),self.Lambda.get())
        out_tomo_path = self._getExtraPath('denoisedEED_'+inp_tomo_path.split('/')[-1])
        args = '{} {} {}'.format(params, inp_tomo_path, out_tomo_path)
        self.runJob('tomoeed', args)
        return out_tomo_path


