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
from pyworkflow import BETA
from pyworkflow.utils import makePath

from tomo3d import Plugin
from tomo.protocols import ProtTomoBase

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, FloatParam

import imod.utils as utils
import os
from tomo.convert import writeTiStack
import tomo.objects as tomoObj
from imod.utils import formatTransformFile


class ProtJjsoftAlignTs(EMProtocol, ProtTomoBase):
    """ Aligning the tilt series using the fiducial positions and motion compensation with tomowarpalign.
    Software from : https://sites.google.com/site/3demimageprocessing/
    Returns the set of aligned tilt series
    """
    _label = 'motion compensated alignment'
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label='Input')
        form.addParam('inputSetOfTiltSeries', PointerParam, important=True,
                      pointerClass='SetOfTiltSeries',
                      label='Input Tilt Series Non Interpolated')

        form.addParam('inputSetOfLandmarkModels', PointerParam, important=True,
                      pointerClass='SetOfLandmarkModels',
                      label='Input Fiducial Models')

        form.addParam('binning', FloatParam,
                      default=1.0,
                      label='Binning',
                      help='Binning to be applied to the interpolated tilt-series. '
                           'Must be a integer bigger than 1')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert every step of the protocol"""
        pre1 = []
        stepId = self._insertFunctionStep(self.convertInputStep)
        pre1.append(stepId)

        pre2 = []
        stepId = self._insertFunctionStep(self.alignTsStep, prerequisites=pre1)
        pre2.append(stepId)

        self._insertFunctionStep(self.computeInterpolatedStackStep, prerequisites=pre2)

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            workingFolder = self._getExtraPath(tsId)
            prefix = os.path.join(workingFolder, tsId)
            makePath(workingFolder)
            tiList = [ti.clone() for ti in ts]
            tiList.sort(key=lambda ti: ti.getTiltAngle())
            tiList.reverse()

            # Creates the st and tlt files
            writeTiStack(tiList,
                         outputStackFn=prefix + '.st',
                         outputTltFn=prefix + '.tlt')

            # Creates the tltxf file
            if ts.getFirstItem().hasTransform():
                formatTransformFile(ts, prefix + '.tltxf')

            # Creates the prexg file as a identity matrix
            self.write_prexg_identity(ts, prefix + '.prexg')

    def alignTsStep(self):
        for fidu in self.inputSetOfLandmarkModels.get():
            tsId = fidu.getTsId()
            prevwFolder = '/'.join(fidu.getFileName().split('/')[:-1])
            workingFolder = self._getExtraPath(tsId)

            aligncom, newstcom = self.get_IMOD_files(prevwFolder, workingFolder, tsId)

            params = ''
            args = '-a {} -n {}'.format(aligncom, newstcom)
            args += params

            self.runJob(Plugin.getTomowarpalignProgram(), args)

    def computeInterpolatedStackStep(self):
        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            extraPrefix = self._getExtraPath(tsId)

            # Naming output tilt series as .mrc
            args = '{}.st {}.mrc'.format(extraPrefix + '/' + tsId, extraPrefix + '/' + tsId)
            self.runJob('cp', args)

            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)
            outputInterpolatedSetOfTiltSeries.append(newTs)

            tltFileName = tsId + ".tlt"
            tltFilePath = os.path.join(self._getExtraPath(tsId), tltFileName)
            tltList = utils.formatAngleList(tltFilePath)

            for index, ti in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(ti, copyId=True)
                newTi.setLocation(index + 1, os.path.join(extraPrefix, '%s.st' % tsId))
                newTi.setTiltAngle(float(tltList[index]))

                if self.binning > 1:
                    newTi.setSamplingRate(ti.getSamplingRate() * int(self.binning.get()))
                newTs.append(newTi)

            if self.binning > 1:
                newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))

            newTs.write()
            outputInterpolatedSetOfTiltSeries.update(newTs)
            outputInterpolatedSetOfTiltSeries.write()

        self._store()

    def createOutputStep(self):
        pass

    # --------------------------- UTILS functions --------------------------------------------
    def getOutputInterpolatedSetOfTiltSeries(self):
        if not hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())

            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputInterpolatedSetOfTiltSeries.setSamplingRate(samplingRate)

            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)

        return self.outputInterpolatedSetOfTiltSeries

    def get_IMOD_files(self, prevwFolder, ts_folder, TsId):
        """Returns the path of the Tilt Series and the angles files"""
        prevprefix = os.path.join(prevwFolder, TsId)
        prefix = os.path.join(ts_folder, TsId)

        self.make_resid_file(prevprefix, prefix)
        aligncom = self.make_aligncom(ts_folder, TsId)
        newstcom = self.make_newstcom(ts_folder, TsId)

        return aligncom, newstcom

    @staticmethod
    def make_resid_file(prevprefix, prefix):
        """Creates the resid file in the correct format from a sfid text file"""
        resid_path = prefix + '.resid'
        sfid_path = prevprefix + '_noGaps.sfid'

        with open(sfid_path)as filex:
            nresids = len(filex.readlines())

        with open(resid_path, 'w') as f:
            f.write('  {} residuals\n'.format(nresids - 1))

            with open(sfid_path) as filex:
                filex.readline()

                for line in filex:
                    line = line.split()
                    f.write('{:10.2f}{:10.2f}{}{:8.2f}{:8.2f}\n'.format(float(line[0]), float(line[1]),
                                                                        str(line[2]).rjust(5),
                                                                        float(line[4]), float(line[5])))
        return resid_path

    @staticmethod
    def write_prexg_identity(ts, prexgPath):
        """Creates the prexg file as a identity matrix"""
        with open(prexgPath, 'w') as f:
            for i in range(len(ts)):
                f.write('1\t0\t0\t1\t0\t0\n')

    @staticmethod
    def make_aligncom(ts_folder, TsId):
        '''Writes an artificial align.com file'''
        aligncomPath = ts_folder + '/align.com'
        pathi = ts_folder + '/'
        with open(aligncomPath, 'w') as f:
            f.write('#ImageSizeXandY    512,512\n\
ImagesAreBinned    1\n\
OutputResidualFile {}.resid\n\
OutputTransformFile {}.tltxf\n\
$xfproduct -StandardInput\n\
InputFile1 {}.prexg\n\
InputFile2 {}.tltxf'.format(*[pathi + TsId] * 4))

            return aligncomPath

    def make_newstcom(self, ts_folder, TsId):
        '''Writes an artifitial newst.com file'''
        newstcomPath = ts_folder + '/newst.com'
        pathi = ts_folder + '/'

        with open(newstcomPath, 'w') as f:
            f.write('$newstack -StandardInput\n\
InputFile	{}.st\n\
OutputFile	{}.ali\n\
TransformFile	{}.xf\n\
TaperAtFill	1,0\n\
AdjustOrigin	\n\
OffsetsInXandY	0.0,0.0\n\
#DistortionField	.idf\n\
ImagesAreBinned	1.0\n\
BinByFactor	1\n\
#GradientFile	{}.maggrad\n\
$if (-e ./savework) ./savework'.format(*[pathi + TsId] * 4))

            return newstcomPath

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        pass

    def _methods(self):
        pass

    def _citations(self):
        return ['Fernandez2018', 'Fernandez2009']

