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
from tomo3d import Plugin
from tomo3d.protocols.protocol_base import ProtBaseTomo3d
from pyworkflow import BETA
from pyworkflow.utils import makePath
from pyworkflow.protocol.params import IntParam, EnumParam, PointerParam, FloatParam, LEVEL_ADVANCED, BooleanParam
from imod.utils import genXfFile

# Motion modelling labels
POLYNOMIAL = 0
SPLINES = 1

# Sample thickness labels
THIN = 0
THICK = 1

# Weighting methods labels
W_NONE = 0
W_RAMP = 1
W_HAMMING = 2
W_SIRT = 3


class ProtJjsoftAlignTomo3dTomogram(ProtBaseTomo3d):
    """ Reconstruct tomograms by aligning the tilt series using the fiducial positions with tomoalign
    and then reconstructs the tomogram with tomorec.
    Software from : https://sites.google.com/site/3demimageprocessing/
    Returns the set of tomograms
    """

    """
        Motion Compensated Reconstruction (ProtJjsoftAlignTomo3dTomogram) — User Manual
            Overview

            The Motion Compensated Reconstruction protocol aligns tilt series using fiducial
            markers and reconstructs tomograms through the tomoalign and tomorec programs.
            Its main objective is to compensate for sample motion and improve tomogram quality
            during cryo-electron tomography reconstruction. By combining fiducial-based alignment
            with weighted reconstruction methods, the protocol generates tomograms with improved
            structural consistency and reduced reconstruction artifacts.

            In biological cryo-ET workflows, this protocol is especially useful for datasets
            affected by beam-induced motion, sample deformation, or acquisition instabilities.
            The use of fiducial landmarks allows more accurate estimation of tilt-series alignment,
            which is critical for preserving structural details in reconstructed tomograms.

            Inputs and General Workflow

            The protocol requires a set of tilt series together with their corresponding fiducial
            landmark models. During preprocessing, tilt angles and fiducial coordinates are converted
            into formats compatible with tomoalign and tomorec. Optional IMOD transformation files
            can also be incorporated to provide initial alignment parameters such as rotations or
            magnification corrections.

            Each tilt series is processed independently. First, the fiducial information is used
            to estimate motion-compensated alignment parameters. Afterwards, the corrected tilt
            series is reconstructed into a tomogram using the selected weighting strategy. Finally,
            the reconstructed tomogram is rotated to restore the correct orientation expected in
            downstream cryo-ET analysis workflows.

            Motion Modeling and Sample Thickness

            The protocol supports two motion modeling strategies: polynomial modeling and spline-based
            interpolation. Polynomial models provide a robust approximation of motion trajectories,
            while splines allow smoother local motion estimation and may better capture complex
            deformations.

            Users can also specify whether the sample is considered thin or thick. Thin samples
            use simpler bivariate models, whereas thick samples require more computationally
            demanding trivariate models to properly account for depth-dependent motion effects.
            Correct selection of sample thickness is important because inaccurate modeling may
            reduce alignment precision and compromise reconstruction quality.

            Reconstruction and Weighting Methods

            After alignment, tomograms are reconstructed using different weighting approaches.
            Ramp and Hamming weighted back-projection methods provide fast reconstruction with
            different levels of noise suppression, while SIRT reconstruction offers iterative
            refinement that can improve contrast and reduce reconstruction artifacts at the cost
            of higher computational time.

            The protocol also allows users to control tomogram dimensions and reconstruction
            boundaries, including width, height, and slice ranges. These options are useful
            for focusing reconstruction on biologically relevant regions while reducing memory
            and computational requirements.

            Outputs and Interpretation

            The protocol generates motion-compensated tomograms for each processed tilt series.
            The reconstructed volumes preserve the spatial information corrected through fiducial
            alignment and can be directly used for visualization, segmentation, subtomogram
            averaging, or downstream structural analysis.

            Because the reconstruction quality strongly depends on fiducial accuracy and motion
            estimation, users should visually inspect the resulting tomograms to verify alignment
            consistency and structural preservation.

            Final Perspective

            Motion-compensated reconstruction is a biologically important step in cryo-electron
            tomography because alignment inaccuracies directly affect the interpretability of
            reconstructed cellular structures. Proper fiducial tracking, appropriate motion
            modeling, and careful selection of reconstruction weighting methods are essential
            for obtaining reliable tomograms suitable for high-quality structural analysis.
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
                      help='Thick samples (thickness > ∼200 nm) require more complex calculations than the thin ones '
                           '(thickness < ∼150 nm).'
                           'Using polynomial motion:\n'
                           '\tthin:  second order bivariate  polynomials\n'
                           '\tthick: second order trivariate polynomials.\n'
                           'Using interpolating splines:\n'
                           '\tthin:  bivariate interpolating splines.\n'
                           '\tthick: trivariate interpolating splines.')
        form.addParam('imodXF', BooleanParam,
                      default=True,
                      label='IMOD Transform file (.xf)',
                      expertLevel=LEVEL_ADVANCED,
                      help='Input IMOD Transform file (.xf) to set the initial alignment parameters (rotation, '
                           'magnification).')

        form.addSection(label='Reconstruction')
        form.addParam('weighting', EnumParam,
                      choices=['None', 'WBP-Ramp', 'WBP-Hamming', 'SIRT'],
                      default=W_HAMMING,
                      label='Weighting Methods',
                      display=EnumParam.DISPLAY_HLIST,
                      help='Methods available for weighting')
        form.addParam('sirtIter', IntParam,
                      default=30,
                      condition='weighting == %i' % W_SIRT,
                      label='SIRT iterations',
                      help='Number of iterations of the SIRT weighting method to be applied')
        self._defineSetShapeParams(form)
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        for ts in self.inputSetOfTiltSeries.get():
            tsId = ts.getTsId()
            fname = ts.getFirstItem().getFileName()
            workingFolder = self.getWorkingDirName(tsId)
            self._insertFunctionStep(self.alignTsStep, tsId, workingFolder)
            self._insertFunctionStep(self.reconstructTomogramStep, tsId, workingFolder, fname)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        for ts, fm in zip(self.inputSetOfTiltSeries.get(), self.inputSetOfLandmarkModels.get()):
            tsId = ts.getTsId()
            workingFolder = self.getWorkingDirName(tsId)
            makePath(workingFolder)
            # Tilt series convert
            if ts.getFirstItem().hasTransform():
                genXfFile(ts, self.getImodXfFile(workingFolder, tsId))
            ts.generateTltFile(self.getAnglesFile(workingFolder, tsId))
            # Fiducials convert
            imodFiducial = self.getImodTxtFiducialsFile(workingFolder, tsId)
            self.parseFiducialFile(fm.getFileName(), imodFiducial)

    def alignTsStep(self, tsId, workingFolder):
        binningFactor = self.binningFactor.get()
        outFile = self.getTomoAlignOutputFile(workingFolder, tsId)
        fiducialTxt = self.getImodTxtFiducialsFile(workingFolder, tsId)
        anglesPath = self.getAnglesFile(workingFolder, tsId)
        transformPath = self.getImodXfFile(workingFolder, tsId)

        params = '-i %s -a %s -o %s ' % (fiducialTxt, anglesPath, outFile)
        if binningFactor > 1:
            params += '-b %.1f ' % binningFactor
        if self.motionModeling.get() == SPLINES:
            params += ' -s '
        params += ' -t %s ' % ('thin' if self.sampleThickness.get() == THIN else 'thick')
        if self.imodXF.get() and exists(transformPath):
            params += ' -I %s ' % transformPath

        self.runJob(Plugin.getTomoAlignProgram(), params)

    def reconstructTomogramStep(self, tsId, workingFolder, tsFileName):
        TsPath, AnglesPath, transformPath = self.getTsFilesMotComp(workingFolder, tsId)
        out_tomo_path = workingFolder + '/tomo_{}.mrc'.format(tsId)
        align_bin = workingFolder + '/alignment_{}.par'.format(tsId)

        params = '-a {} -i {} -o {}'.format(align_bin, tsFileName, out_tomo_path)
        if self.weighting.get() == W_RAMP:
            params += ' -w ramp'
        elif self.weighting.get() == W_HAMMING:
            params += ' -w hamming'
        elif self.weighting.get() == W_SIRT:
            params += ' -w sirt -l %i' % self.sirtIter.get()
        else:
            params += ' -w none'

        if self.setShape.get() == 0:
            if self.width.get() != 0:
                params += ' -x {}'.format(self.width.get())
            if self.finSlice.get() != 0:
                params += ' -Y {},{}'.format(self.iniSlice.get(), self.finSlice.get())
            if self.height.get() != 0:
                params += ' -z {}'.format(self.height.get())

        if self.imodXF.get() and exists(transformPath):
            params += ' -S ' + transformPath

        params += ' -t %i ' % self.numberOfThreads.get()

        self.runJob(Plugin.getTomoRecProgram(), params)
        out_tomo_rx_path = self.rotXTomo(tsId)
        self.outputFiles.append(out_tomo_rx_path)

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        pass

    def _methods(self):
        pass

    # --------------------------- UTILS functions --------------------------------------------
    @staticmethod
    def getTsFilesMotComp(tsFolder, tSId):
        """Returns the path of the Tilt Serie and the angles files"""
        prefix = join(tsFolder, tSId)
        TsPath = prefix + '.st'
        AnglesPath = prefix + '.rawtlt'
        transformPath = prefix + '.xf'

        return TsPath, AnglesPath, transformPath

    @staticmethod
    def parseFiducialFile(scipionFid, outFid):
        """Converts the Scipion fid format to JJ format needed"""
        #TODO: parse with csvreader and substract 1 to Z column
        with open(outFid, 'w') as f:
            with open(scipionFid) as filex:
                filex.readline()
                for line in filex:
                    p = line.split('\t')
                    f.write('1\t{}\t{}\t{}\t{}\n'.format(p[3], float(p[0]), float(p[1]), float(p[2])))

    @staticmethod
    def getPathAndBaseName(workingFolder, tsId):
        return join(workingFolder, tsId)

    @staticmethod
    def getAnglesFile(workingFolder, tsId):
        return ProtJjsoftAlignTomo3dTomogram.getPathAndBaseName(workingFolder, tsId) + '.tlt'

    @staticmethod
    def getImodXfFile(workingFolder, tsId):
        return ProtJjsoftAlignTomo3dTomogram.getPathAndBaseName(workingFolder, tsId) + '.xf'

    @staticmethod
    def getImodTxtFiducialsFile(workingFolder, tsId):
        return join(workingFolder, 'imod_%s.fid.txt' % tsId)

    @staticmethod
    def getTomoAlignOutputFile(workingFolder, tsId):
        return join(workingFolder, 'alignment_%s.par' % tsId)
