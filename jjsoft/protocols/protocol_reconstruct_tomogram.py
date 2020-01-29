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
from pyworkflow.protocol.params import IntParam, EnumParam, PointerParam, FloatParam

from tomo.objects import Tomogram
import os
import pyworkflow as pw
from tomo.convert import writeTiStack

class JjsoftReconstructTomogram(EMProtocol, ProtTomoBase):
    """ Reconstruct tomograms from aligned tilt series using TOMO3D from
    Software from: https://sites.google.com/site/3demimageprocessing/
    Returns the set of tomograms
    """
    _label = 'reconstruct tomogram'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label='Input')
        form.addParam('inputSetOfTiltSeries', PointerParam, important=True,
                      pointerClass='SetOfTiltSeries',
                      label='Input Tilt Series')
        form.addParam('method', EnumParam,
                      choices=['WBP (Fast)','SIRT (Slow)'], default=0,
                      label='Reconstruction method',
                      help='Reconstrution method to use')
        form.addParam('nIterations', IntParam, default=30,
                      condition='method==1',
                      label='Number of Iterations (SIRT)',
                      help='Number of Iterations used in the SIRT method')
        form.addParam('Hamming', FloatParam, default=0.0,
                      label='Hamming filter frequency',
                      help='Frequency for the Hamming atenuation filter [0,0.5]. \n0 always uses the filter, 0.5 turns it off')

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
        group.addParam('height', IntParam,
                       default=0,
                       label='Thickness',
                       help='Height of the reconstructed tomogram (Default: width of the tomogram)')
        group.addParam('iniSlice', IntParam,
                       default=0,
                       label='Initial slice',
                       help='Initial slice (of range) to include')
        group.addParam('finSlice', IntParam,
                       default=0,
                       label='Final slice',
                       help='Final slice (of range) to include (Maximum must be the size of tilt series)')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert every step of the protocol"""
        pre1 = []
        stepId = self._insertFunctionStep('convertInputStep')
        pre1.append(stepId)

        pre = []
        self.outputFiles = []
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            stepId = self._insertFunctionStep('reconstructTomogramStep', tsId, workingFolder, prerequisites=pre1)
            pre.append(stepId)

        self._insertFunctionStep('createOutputStep', prerequisites=pre)

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
                         outputTltFn=prefix + '.rawtlt')

    def reconstructTomogramStep(self, tsId, workingFolder):
        # We start preparing writing those elements we're using as input to keep them untouched
        TsPath, AnglesPath = self.get_Ts_files(workingFolder, tsId)
        out_tomo_path = workingFolder + '/tomo_{}.mrc'.format(tsId)
        params=''
        if self.method==1:
            params+=' -S -l '+str(self.nIterations)
        if self.setShape.get() == 0:
            if self.width.get() != 0:
                params += ' -x {}'.format(self.width.get())
            if self.finSlice.get() != 0:
                params += ' -y {},{}'.format(self.iniSlice.get(), self.finSlice.get())
            if self.height.get() != 0:
                params += ' -z {}'.format(self.height.get())

        args = '-i {} -a {} -o {} -t {}'.format(TsPath, AnglesPath, out_tomo_path, self.numberOfThreads)
        args+=params
        self.runJob('tomo3d', args)
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
        AnglesPath = prefix + '.rawtlt'
        return TsPath, AnglesPath



