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
from os.path import join

import mrcfile
import numpy as np
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.utils import makePath

from jjsoft import Plugin
from tomo.protocols import ProtTomoBase

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import IntParam, EnumParam, PointerParam, FloatParam

from tomo.objects import Tomogram
import os


class ProtBaseReconstruct(EMProtocol, ProtTomoBase):
    """ Reconstruct tomograms from aligned tilt series using TOMO3D from
    Software from: https://sites.google.com/site/3demimageprocessing/
    Returns the set of tomograms
    """

    # --------------------------- DEFINE param functions --------------------------------------------
    # def __init__(self, **kwargs):
    #     super().__init__(**kwargs)
        # workingFolder = None

    def _defineParams(self, form):
        pass

    @staticmethod
    def _defineInputParams(form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label='Input')
        form.addParam('inputSetOfTiltSeries', PointerParam, important=True,
                      pointerClass='SetOfTiltSeries',
                      label='Non Interpolated Tilt Series')

    @staticmethod
    def _defineReconstructParams(form):
        form.addParam('inputSetOfTiltSeries', PointerParam, important=True,
                      pointerClass='SetOfTiltSeries',
                      label='Input Tilt Series')
        form.addParam('method', EnumParam,
                      choices=['WBP (Fast)', 'SIRT (Slow)'], default=0,
                      label='Reconstruction method',
                      help='Reconstrution method to use')
        form.addParam('nIterations', IntParam, default=30,
                      condition='method==1',
                      label='Number of Iterations (SIRT)',
                      help='Number of Iterations used in the SIRT method')
        form.addParam('Hamming', FloatParam, default=0.0,
                      label='Hamming filter frequency',
                      help='Frequency for the Hamming atenuation filter [0,0.5]. \n0 always uses the filter, '
                           '0.5 turns it off')
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
        pass

    # --------------------------- STEPS functions --------------------------------------------
    def reconstructTomogramStep(self, tsId, workingFolder):
        # We start preparing writing those elements we're using as input to keep them untouched
        TsPath, AnglesPath = self.getTsFiles(workingFolder, tsId)
        out_tomo_path = workingFolder + '/tomo_{}.mrc'.format(tsId)
        params = ''
        if self.method == 1:
            params += ' -S -l '+str(self.nIterations)
        if self.setShape.get() == 0:
            if self.width.get() != 0:
                params += ' -x {}'.format(self.width.get())
            if self.finSlice.get() != 0:
                params += ' -y {},{}'.format(self.iniSlice.get(), self.finSlice.get())
            if self.height.get() != 0:
                params += ' -z {}'.format(self.height.get())

        args = '-i {} -a {} -o {} -t {}'.format(TsPath, AnglesPath, out_tomo_path, self.numberOfThreads)
        args += params
        self.runJob(Plugin.getTomoRecProgram(), args)
        out_tomo_rx_path = self.rotXTomo(tsId)
        self.outputFiles.append(out_tomo_rx_path)

    def rotXTomo(self, tsId):
        """Result of the reconstruction must be rotated 90 degrees around the X
        axis to recover the original orientation (due to jjsoft design)"""
        outPath = self._getExtraPath(tsId)
        makePath(outPath)
        inTomoFile = join(self._getTmpPath(tsId), 'tomo_%s.mrc' % tsId)
        outTomoFile = join(outPath, 'tomo_%s.mrc' % tsId)

        with mrcfile.open(inTomoFile, mode='r', permissive=True) as mrc:
            rotData = np.rot90(mrc.data)

        with mrcfile.open(outTomoFile, mode='w+') as mrc:
            mrc.set_data(rotData)

        return outTomoFile

    def createOutputStep(self):
        outputTomos = self._createSetOfTomograms()
        outputTomos.copyInfo(self.inputSetOfTiltSeries.get())

        for tomoPath, ts in zip(self.outputFiles, self.inputSetOfTiltSeries.get()):
            tomo = Tomogram()
            tomo.setLocation(tomoPath)
            tomo.setSamplingRate(ts.getSamplingRate())
            tomo.setOrigin()
            tomo.setTsId(ts.getTsId())
            outputTomos.append(tomo)

        self._defineOutputs(outputTomograms=outputTomos)
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
    @staticmethod
    def getTsFiles(tsFolder, tsId):
        """Returns the path of the Tilt Serie and the angles files"""
        prefix = join(tsFolder, tsId)
        TsPath = prefix + '.st'
        AnglesPath = prefix + '.rawtlt'
        return TsPath, AnglesPath



