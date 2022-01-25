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
from os.path import join, exists

import mrcfile
import numpy as np

from jjsoft import Plugin
from jjsoft.protocols.protocol_base_reconstruct import ProtBaseReconstruct
from pyworkflow import BETA
from pyworkflow.utils import makePath
from pyworkflow.protocol.params import IntParam, EnumParam, PointerParam, FloatParam, LEVEL_ADVANCED, BooleanParam
from tomo.objects import Tomogram
from imod.utils import formatTransformFile

# Motion modelling labels
POLYNOMIAL = 0
SPLINES = 1

# Sample thickness labels
THIN = 0
THICK = 1

# Weighting methods labels
# ['None', 'WBP-Ramp', 'WBP-Hamming', 'SIRT']
W_NONE = 0
W_RAMP = 1
W_HAMMING = 2
W_SIRT = 3


class ProtJjsoftAlignReconstructTomogram(ProtBaseReconstruct):
    """ Reconstruct tomograms by aligning the tilt series using the fiducial positions with tomoalign
    and then reconstructs the tomogram with tomorec.
    Software from : https://sites.google.com/site/3demimageprocessing/
    Returns the set of tomograms
    """
    _label = 'motion compensated reconstruction'
    _devStatus = BETA

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        # First we customize the inputParticles param to fit our needs in this protocol
        self._defineInputParams(form)
        form.addParam('inputSetOfLandmarkModels', PointerParam, important=True,
                      pointerClass='SetOfLandmarkModels',
                      label='Input Fiducial Models')

        form.addSection(label='Alignment')
        form.addParam('binningFactor', FloatParam,
                      default=1.0,
                      label='Binning factor',
                      help='Binning to be applied to the interpolated tilt-series.')
        form.addParam('motionModeling', EnumParam,
                      choices=['Polynomial', 'Splines'],
                      default=POLYNOMIAL,
                      label='Motion modeling',
                      display=EnumParam.DISPLAY_HLIST,
                      help='Motion modelling by polynomials or by interpolating splines.')
        form.addParam('sampleThickness', EnumParam,
                      choices=['Thin', 'Thick'],
                      default=THIN,
                      label='Sample Thickness',
                      display=EnumParam.DISPLAY_HLIST,
                      help='Thick samples require more complex calculations than thinner.'
                           'Using polynomial motion:\n'
                           '\tthin:  second order bivariate  polynomials\n'
                           '\tthick: second order trivariate polynomials.\n'
                           'Using interpolating splines:\n'
                           '\tthin:  bivariate interpolating splines.\n'
                           '\tthick: trivariate interpolating splines.')
        form.addParam('imodXF', BooleanParam,
                      default=True,
                      label='IMOD Transform file (.xf)',
                      display=EnumParam.DISPLAY_HLIST,
                      expertLevel=LEVEL_ADVANCED,
                      help='Input IMOD Transform file (.xf) to set the initial alignment parameters (rotation, '
                           'magnification).')

        form.addSection(label='Reconstruction')
        form.addParam('weighting', BooleanParam,
                      default=True,
                      label='Apply weighting?',
                      help='Methods available for weighting')
        form.addParam('sirtIter', IntParam,
                      default=30, condition='weighting==3',
                      label='SIRT iterations',
                      help='Number of iterations of the SIRT weighting method to be applied')

        self._defineReconstructParams(form)
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        for ts, fm in zip(self.inputSetOfTiltSeries.get(), self.inputSetOfLandmarkModels.get()):
            tsId = ts.getTsId()
            self.workingFolder = self._getExtraPath(tsId, fm.getFileName())
            makePath(self.workingFolder)
            self._insertFunctionStep(self.convertInputStep, tsId)
            self._insertFunctionStep(self.alignTsStep, tsId)
            self._insertFunctionStep(self.reconstructTomogramStep, tsId)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self, tsId, scipionFiducial):
        # Tilt series convert
        prefix = join(self.workingFolder, tsId)
        if self.ts.getFirstItem().hasTransform():
            formatTransformFile(self.ts, prefix + '.xf')

        # Fiducial model convert
        imodFiducial = join(self.workingFolder, 'imod_%s.fid.txt' % tsId)
        self.parse_fid(scipionFiducial, imodFiducial)

    def alignTsStep(self, tsId):
        TsPath, AnglesPath, transformPath = self.getTsFiles(self.workingFolder, tsId)
        fiducial_text = join(self.workingFolder, 'imod_%s.fid.txt' % tsId)
        out_bin = join(self.workingFolder, 'alignment_%s.bin' % tsId)
        binningFactor = self.binningFactor.get()

        params = '-i %s -a %s -o %s' % (fiducial_text, AnglesPath, out_bin)
        if binningFactor > 1:
            params += '-b %.1f ' % binningFactor
        if self.motionModeling.get() == SPLINES:
            params += '-s '
        params += '-t %s ' % ('thin' if self.sampleThickness.get() == THIN else 'thick')
        if self.imodXF.get() and exists(transformPath):
            params += ' -I %s ' % transformPath

        self.runJob(Plugin.getTomoAlignProgram(), params)

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

    # --------------------------- UTILS functions --------------------------------------------
    @staticmethod
    def getTsFiles(tsFolder, tSId):
        """Returns the path of the Tilt Serie and the angles files"""
        prefix = join(tsFolder, tSId)
        TsPath = prefix + '.st'
        AnglesPath = prefix + '.tlt'
        transformPath = prefix + '.xf'

        return TsPath, AnglesPath, transformPath

    @staticmethod
    def parse_fid(scip_fid,out_fid):
        """Converts the scipion fid format to JJ format needed"""
        with open(out_fid, 'w') as f:
            with open(scip_fid) as filex:
                filex.readline()
                for line in filex:
                    p = line.split('\t')
                    f.write('1\t{}\t{}\t{}\t{}\n'.format(p[3], float(p[0]), float(p[1]), float(p[2])))




