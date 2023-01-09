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
from pyworkflow.protocol.params import IntParam, EnumParam, FloatParam

import os

# Reconstruction methods
WBP = 0
SIRT = 1


class ProtJjsoftReconstructTomogram(ProtBaseReconstruct):
    """ Reconstruct tomograms using TOMO3D from
    Software from: https://sites.google.com/site/3demimageprocessing/
    Returns the set of tomograms
    """
    _label = 'reconstruct tomogram'
    _devStatus = BETA

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineInputParams(form)
        form.addParam('method', EnumParam,
                      choices=['WBP (Fast)', 'SIRT (Slow)'], default=WBP,
                      label='Reconstruction method',
                      help='Reconstrution method to use')
        form.addParam('nIterations', IntParam, default=30,
                      condition='method==%i' % SIRT,
                      label='Number of Iterations (SIRT)',
                      help='Number of Iterations used in the SIRT method')
        form.addParam('Hamming', FloatParam, default=0.0,
                      label='Hamming filter frequency',
                      help='Frequency for the Hamming atenuation filter [0,0.5]. \n0 always uses the filter, '
                           '0.5 turns it off')
        self._defineSetShapeParams(form)
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

            ts.applyTransform(outputStackFn)
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


