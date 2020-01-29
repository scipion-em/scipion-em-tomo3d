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
from tomo.protocols import ProtTomoBase

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import IntParam, EnumParam, PointerParam, FloatParam, LEVEL_ADVANCED

from tomo.objects import Tomogram
import os
import pyworkflow as pw
from tomo.convert import writeTiStack

class JjsoftAlignReconstructTomogram(EMProtocol, ProtTomoBase):
    """ Reconstruct tomograms by aligning the tilt series using the fiducial positions with tomoalign
    and then reconstructs the tomogram with tomorec. It returns the set of tomograms
    """
    _label = 'align and reconstruct tomogram'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label='Input')
        form.addParam('inputSetOfTiltSeries', PointerParam, important=True,
                      pointerClass='SetOfTiltSeries',
                      label='Input Tilt Series')
        form.addParam('inputSetOfLandmarkModels', PointerParam, important=True,
                      pointerClass='SetOfLandmarkModels',
                      label='Input Fiducial Models')

        form.addSection(label='Alignment')
        form.addParam('useBinning', EnumParam,
                      choices=['Yes', 'No'],
                      default=0,
                      label='Use binning',
                      display=EnumParam.DISPLAY_HLIST,
                      help='Use of binning factor to be applied internally to the input fiducial file')
        form.addParam('binningFactor', FloatParam,
                       default=1.0, condition = 'useBinning==0',
                       label='Binning factor',
                       help='Binning to be applied to the interpolated tilt-series. '
                            'Must be a integer bigger than 1')
        form.addParam('motionModeling', EnumParam,
                      choices=['Polynomial', 'Interpolation'],
                      default=0,
                      label='Motion modeling',
                      display=EnumParam.DISPLAY_HLIST,
                      help='Modeling of motion by polynomials or by interpolation')
        form.addParam('sampleThickness', EnumParam,
                      choices=['Thin', 'Thick'],
                      default=0,
                      label='Sample Thickness',
                      display=EnumParam.DISPLAY_HLIST,
                      help='The Z direction will only be taken into account in thick samples')

        form.addSection(label='Reconstruction')
        form.addParam('weighting', EnumParam,
                      choices=['None', 'WBP-Ramp','WBP-Hamming','SIRT'],
                      default=2,
                      label='Weighting Methods',
                      display=EnumParam.DISPLAY_HLIST,
                      help='Methods available for weighting')
        form.addParam('sirtIter', IntParam,
                      default=30, condition='weighting==3',
                      label='SIRT iterations',
                      help='Number of iterations of the SIRT weighting method to be applied')
        form.addParam('imodXF', EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Use IMOD xf',
                      display=EnumParam.DISPLAY_HLIST,
                      help='Use transformation file xf from IMOD to apply in the original stack',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('setShape', EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Set manual tomogram shape',
                      display=EnumParam.DISPLAY_HLIST,
                      help='By deafault the shape of the tomogram is defined by the tilt series shape')

        group = form.addGroup('Tomogram shape', condition='setShape==0')
        group.addParam('width', IntParam,
                       default=0,
                       label='Width',
                       help='Focus the tomogram in a region of the tilt series')
        group.addParam('nSlices', IntParam,
                       default=0,
                       label='Number of slices',
                       help='The number of the slices (along the tilt axis) of the reconstructed tomogram')
        group.addParam('height', IntParam,
                       default=0,
                       label='Thickness',
                       help='Height of the reconstructed tomogram (Default: width of the tomogram)')


        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert every step of the protocol"""
        pre1 = []
        stepId = self._insertFunctionStep('convertInputStep')
        pre1.append(stepId)

        pre2,pre3,pre4 = [], [], []
        self.outputFiles = []

        stepId = self._insertFunctionStep('convertFiducialTextStep', prerequisites=pre1)
        pre2.append(stepId)
        stepId = self._insertFunctionStep('alignTsStep', prerequisites=pre2)
        pre3.append(stepId)
        stepId = self._insertFunctionStep('reconstructTomogramStep', prerequisites=pre3)
        pre4.append(stepId)

        self._insertFunctionStep('createOutputStep', prerequisites=pre4)

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            prefix = os.path.join(workingFolder, tsId)
            pw.utils.makePath(workingFolder)
            tiList = [ti.clone() for ti in ts]
            tiList.sort(key=lambda ti: ti.getTiltAngle())
            tiList.reverse()
            writeTiStack(tiList,
                         outputStackFn=prefix + '.st',
                         outputTltFn=prefix + '.tlt')

    def convertFiducialTextStep(self):
        for fidu in self.inputSetOfLandmarkModels.get():
            tsId = fidu.getTsId()
            workingFolder = self._getExtraPath(tsId)

            scip_fiducial = fidu.getFileName()
            imod_fiducial = workingFolder + '/imod_{}.fid.txt'.format(tsId)

            self.parse_fid(scip_fiducial, imod_fiducial)

    def alignTsStep(self):
        for fidu in self.inputSetOfLandmarkModels.get():
            tsId = fidu.getTsId()
            workingFolder = self._getExtraPath(tsId)
            TsPath, AnglesPath = self.get_Ts_files(workingFolder, tsId)
            fiducial_text = workingFolder + '/imod_{}.fid.txt'.format(tsId)
            out_bin = workingFolder + '/alignment_{}.bin'.format(tsId)

            params = ''
            if self.useBinning.get() == 0:
                params += ' -b {}'.format(self.binningFactor.get())
            if self.motionModeling.get() == 1:
                params += ' -s'
            if self.sampleThickness.get() == 0:
                params += ' -t thin'
            else:
                params += ' -t thick'

            args = '-i {} -a {} -o {}'.format(fiducial_text, AnglesPath, out_bin)
            args += params
            self.runJob('tomoalign', args)

    def reconstructTomogramStep(self):
        for fidu in self.inputSetOfLandmarkModels.get():
            tsId = fidu.getTsId()
            workingFolder = self._getExtraPath(tsId)
            TsPath, AnglesPath = self.get_Ts_files(workingFolder, tsId)
            out_tomo_path = workingFolder + '/tomo_{}.mrc'.format(tsId)
            align_bin = workingFolder + '/alignment_{}.bin'.format(tsId)

            params = ''
            if self.weighting.get() == 1:
                params += ' -w ramp'
            elif self.weighting.get() == 2:
                params += ' -w hamming'
            elif self.weighting.get() == 3:
                params += ' -w sirt -l {}'.format(self.sirtIter.get())

            if self.setShape.get() == 0:
                if self.width.get() != 0:
                    params += ' -x {}'.format(self.width.get())
                if self.nSlices.get() != 0:
                    params += ' -y {}'.format(self.nSlices.get())
                if self.height.get() != 0:
                    params += ' -z {}'.format(self.height.get())

            if self.imodXF.get() == 0:
                xf_file = fidu.getFileName()
                xf_file = '/'.join(xf_file.split('/')[:-1]) + '/{}local.xf'.format(tsId)
                params += ' -S '+xf_file

            params += ' -t {}'.format(self.numberOfThreads)

            args = '-a {} -i {} -o {}'.format(align_bin, TsPath, out_tomo_path)
            args += params
            self.runJob('tomorec', args)
            self.outputFiles.append(out_tomo_path)

    def createOutputStep(self):
        outputTomos = self._createSetOfTomograms()
        outputTomos.copyInfo(self.inputSetOfTiltSeries.get())
        #outputTomos.setSamplingRate(self.inputSetOfTiltSeries.getSamplingRate())

        for i, inp_ts in enumerate(self.inputSetOfTiltSeries.get()):
            tomo_path = self.outputFiles[i]
            tomo = Tomogram()
            tomo.setLocation(tomo_path)
            tomo.setSamplingRate(inp_ts.getSamplingRate())
            # tomo.setAcquisition(inp_ts.getAcquisition())
            outputTomos.append(tomo)


        self._defineOutputs(outputTomograms=outputTomos)
        self.outputTomograms=outputTomos
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
    def get_Ts_files(self,ts_folder,TsId):
        '''Returns the path of the Tilt Serie and the angles files'''
        prefix = os.path.join(ts_folder, TsId)
        TsPath = prefix + '.st'
        AnglesPath = prefix + '.tlt'
        return TsPath, AnglesPath

    def parse_fid(self,scip_fid,out_fid):
        '''Converts the scipion fid format to JJ format needed'''
        with open(out_fid, 'w') as f:
            with open(scip_fid) as filex:
                filex.readline()
                for line in filex:
                    p = line.split('\t')
                    f.write('1\t{}\t{}\t{}\t{}\n'.format(p[3],float(p[0]),float(p[1]),float(p[2])))



