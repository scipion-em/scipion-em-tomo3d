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
from os.path import join


from tomo3d.protocols.protocol_base_reconstruct import ProtBaseReconstruct
from pyworkflow import BETA
from pyworkflow.utils import makePath

from tomo3d import Plugin
from pyworkflow.protocol.params import IntParam, EnumParam, FloatParam, LEVEL_ADVANCED

import os

# Reconstruction methods
WBP = 0
SIRT = 1


class ProtJjsoftReconstructTomogram(ProtBaseReconstruct):
    """ Reconstruct tomograms using TOMO3D from
    Software from: https://sites.google.com/site/3demimageprocessing/

    Tomo3D implements a multithreaded vectorized approach to tomographic reconstruction
    that takes full advantage of the power in modern multicore computers. Full resolution
    tomograms are generated at high speed on standard computers with no special system
    requirements. Tomo3D has the most common reconstruction methods implemented,
    namely WBP and SIRT. It proves to be competitive with current GPU solutions in terms
    of processing time, in the order of a few seconds with WBP or minutes with SIRT.
    """
    _label = 'reconstruct tomogram'
    _devStatus = BETA

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineInputParams(form)
        form.addParam('method', EnumParam,
                      choices=['WBP (Fast)', 'SIRT (Slow)'], default=WBP,
                      label='Reconstruction method',
                      help='Reconstrution method to use:\n'
                           '_WBP_: Also called filtered back projection provides an exact'
                           'solution for the reconstruction problem. Unfortunately, the '
                           'contrast of the reconstructed tomogram is poor. This solution should'
                           'be chosen for subtomogram averaging approaches\n'
                           '_SIRT_: It is an algebraic iterative algorithm that aims to maximize'
                           'the compabitity between the images and the reconstructed tomogram. This'
                           'solution is more suitable for celullar environment, but not for '
                           'subtomogram averaging')
        form.addParam('nIterations', IntParam, default=30,
                      condition='method==%i' % SIRT,
                      label='Number of Iterations (SIRT)',
                      help='It sets the number of iterations for SIRT. By default,'
                           ' a number of 30 iterations is used.')
        self._defineSetShapeParams(form)
        form.addParam('Hamming', FloatParam, default=0.0,
                      condition='method==%i' % WBP,
                      expertLevel=LEVEL_ADVANCED,
                      label='Hamming filter frequency',
                      help=' By default, in WBP the program applies a Ramp filter together '
                           ' with a Hamming filter. The latter allows attenuation of the '
                           ' high frequency components, which are typically very noisy. '
                           ' This parameter allows setting up the starting frequency from '
                           ' which the Hamming filter is to be applied. The values '
                           ' must range in [0,0.5], where 0.5 represents the Nyquist frequency.'
                           ' By default, the program uses a starting frequency of 0, i.e. a '
                           ' Hamming is applied over the whole frequency range. A value of 0.5 '
                           ' would turn off the Hamming filter. Values in-between would '
                           ' preserve low-frequency components and modulate the contribution '
                           ' of the other frequency components')
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert every step of the protocol"""
        self._insertFunctionStep(self.convertInputStep)

        for ts in self.inputSetOfTiltSeries.get():

            tsId = ts.getTsId()
            workingFolder = self.getWorkingDirName(tsId)
            self._insertFunctionStep(self.reconstructTomogramStep, tsId, workingFolder)

        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        for ts in self.inputSetOfTiltSeries.get():

            excludedViewsCount = len(ts.getExcludedViewsIndex())

            if excludedViewsCount > 0:
                msg = "Invalid input: Tomo 3d does not work with excluded views and Tilt series %s has %d " \
                      "Try to remove them with imod exclude views protocol. Stopping now." % (ts.getTsId(), excludedViewsCount)
                self.warning(msg)
                raise AttributeError(msg)

            tsId = ts.getTsId()
            workingFolder = self.getWorkingDirName(tsId)
            prefix = os.path.join(workingFolder, tsId)
            makePath(workingFolder)

            outputStackFn = prefix + '.st'
            outputTltFn = prefix + '.rawtlt'

            rotationAngle = ts.getAcquisition().getTiltAxisAngle()

            # Check if rotation angle is greater than 45º. If so,
            # swap x and y dimensions to adapt output image sizes to
            # the final sample disposition.
            swapTrue = False
            if 45 < abs(rotationAngle) < 135:
                swapTrue = True
            ts.applyTransform(outputStackFn, swapTrue)
            ts.generateTltFile(outputTltFn)

    def reconstructTomogramStep(self, tsId, workingFolder):
        # We start preparing writing those elements we're using as input to keep them untouched
        TsPath, AnglesPath = self.getTsFiles(workingFolder, tsId)
        outTomoPath = join(workingFolder, '%s.mrc' % tsId)
        params = ''
        if self.method.get() == SIRT:
            params += ' -S -l %i ' % self.nIterations.get()
        if self.setShape.get():
            if self.width.get() != 0:
                params += ' -x {}'.format(self.width.get())
            if self.finSlice.get() != 0:
                params += ' -y {},{}'.format(self.iniSlice.get(), self.finSlice.get())
            if self.height.get() != 0:
                params += ' -z {}'.format(self.height.get())

        args = '-i {} -a {} -o {} -t {}'.format(TsPath, AnglesPath, outTomoPath, self.numberOfThreads)
        args += params
        self.runJob(Plugin.getTomo3dProgram(), args)
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

    def _citations(self):
        return ['Fernandez2010_tomo3d', 'Fernandez2015_tomo3d']

# --------------------------- UTILS functions --------------------------------------------


