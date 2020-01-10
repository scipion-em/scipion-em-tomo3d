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

    DENOISE_AND = 0
    DENOISE_BF = 1
    DENOISE_EED = 2

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
                      choices=['Anistropic Non-linear Diffusion','BFlow', 'Edge Enhancing Diffusion'], default=0,
                      label='Denoising method',
                      help='Denoising method to use.')
        form.addSection(label='Parameters')
        form.addParam('SigmaGaussian', FloatParam, default=0.5,
                      label='Sigma Gaussian Filter',
                      help='Sigma for initial gaussian filtering.')
        form.addParam('nIter', IntParam, default=40,
                      label='Number of Iterations',
                      help='Number of Iterations of denoising.')
        #
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
            print('Denoising by Anistropic Non-linear Diffusion')
            out_tomo_path = self.call_AND(inp_tomo_path, iter)

        elif  self.method.get() == 1:
            print('Denoising by BFlow')
            # call BFlow
            out_tomo_path = self.call_BFlow(inp_tomo_path)

        elif self.method.get() == 2:
            print('Denoising by Edge Enhancing Diffusion')
            out_tomo_path = self.call_EED(inp_tomo_path)

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
        return ['Fernandaz2018','Fernandez2009','Fernandez2003']

    # --------------------------- UTILS functions --------------------------------------------
    def write_AND_file(self, csh_name, input_tomo_path, output_tomo_path):
        '''Writes the AND.csh file'''
        f=open(csh_name,'w')
        f.write(
            'tomoand << eof\n'
            '{}\t # Input file\n'
            '{}\t # Output file\n'
            '{} N\t # Number of Iterations and [Use of Stopping Criterion]\n'
            '-1\t # C Constant for CED\n'
            '-1\t # K Constant for EED\n'
            '0.5\t # CED/EED Balance Parameter\n'
            '0.5\t # Proportion of CED along 2nd eigenvector\n'
            '0.5\t # Proportion of Smoothing based on Grey level\n'
            '45 140 45\t # Coordinates of the Noise Area ( X Y Z )\n'
            '{} 2\t # Initial sigma and sigma for averaging struc.Tensor\n'
            '{} {}\t # ht and [N. Threads]\n'
            '0\t # [Diffusion Mode: Standard Hybrid (0); Averaged EED(1|2)]\n'
            'eof'.format(input_tomo_path,output_tomo_path,
                        self.nIter.get(),self.SigmaGaussian.get(),
                        self.TimeStep.get(),self.numberOfThreads)
        )
        f.close()
        #return params['output_tomo_path']

    def create_AND_file(self, input_tomo_path, iter):
        '''Create the csh file that stores the parameters'''
        output_csh = self._getExtraPath('tomoand_{}.csh'.format(iter))
        output_tomo_path = self._getExtraPath('denoisedEED_'+input_tomo_path.split('/')[-1])

        self.write_AND_file(output_csh, input_tomo_path, output_tomo_path)

        args=''
        call='chmod 777 {}'.format(output_csh)
        self.runJob(call, args)
        return output_csh, output_tomo_path

    def call_AND(self, inp_tomo_path, iter):
        '''Denoises de tomogram using the AND method'''
        and_csh, out_tomo_path = self.create_AND_file(inp_tomo_path, iter)
        args=''
        self.runJob(and_csh,args)
        return out_tomo_path

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
        params = '-g {} -i {} -s {}'.format(self.SigmaGaussian.get(), self.nIter.get(), self.TimeStep.get())
        out_tomo_path = self._getExtraPath('denoisedEED_'+inp_tomo_path.split('/')[-1])
        args = '{} {} {}'.format(params, inp_tomo_path, out_tomo_path)
        self.runJob('tomoeed', args)
        return out_tomo_path


