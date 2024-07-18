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
from enum import Enum
from os.path import join

import mrcfile
import numpy as np

from pyworkflow.object import Set
from pyworkflow.utils import makePath

from tomo.protocols import ProtTomoBase

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import IntParam, EnumParam, PointerParam, BooleanParam, LEVEL_ADVANCED

from tomo.objects import Tomogram, SetOfTomograms


class outputTomoRecObjects(Enum):
    tomograms = SetOfTomograms


class ProtBaseReconstruct(EMProtocol, ProtTomoBase):
    """ Reconstruct tomograms from aligned tilt series using TOMO3D from
    Software from: https://sites.google.com/site/3demimageprocessing/
    Returns the set of tomograms
    """
    _OUTNAME = outputTomoRecObjects.tomograms.name
    _possibleOutputs = {_OUTNAME: SetOfTomograms}

    # --------------------------- DEFINE param functions --------------------------------------------
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None

    def _defineParams(self, form):
        pass

    @staticmethod
    def _defineInputParams(form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label='Input')
        form.addParam('inputSetOfTiltSeries', PointerParam, important=True,
                      pointerClass='SetOfTiltSeries',
                      label='Tilt series',
                      help='Tilt Series to reconstruct the tomograms. Ideally these tilt series'
                            'should contain alignment information.')

    @staticmethod
    def _defineSetShapeParams(form):
        form.addParam('setShape', BooleanParam,
                      default=True,
                      label='Set manual tomogram shape',
                      display=BooleanParam,
                      help='By default the shape of the tomogram is defined by the tilt series shape,'
                           'but the user can manually set it here')

        group = form.addGroup('Tomogram dimensions', condition='setShape')

        group.addParam('height', IntParam,
                       default=0,
                       label='Thickness (px)',
                       help='By default, the height (i.e. thickness) of the reconstructed tomogram '
                            ' is equal to the X dimension of the images in the tilt-series. '
                            ' This option allows the user to properly set the height of the '
                            ' tomogram to a value that better fits the thickness of the '
                            ' specimen under study. This also allows reduction of the processing '
                            ' time. The thickness value has to be an even number. This is '
                            ' necessary to exploit projection symmetry (see Agulleiro et al. '
                            ' J.Struct.Biol. 170:570–575, 2010).')

        group.addParam('width', IntParam,
                       expertLevel=LEVEL_ADVANCED,
                       default=0,
                       label='Width (px)',
                       help='By default, the width of the reconstructed tomogram is equal to '
                            ' the X dimension of the images in the tilt-series. This option '
                            ' allows the user to set the width of the tomogram to a smaller '
                            ' value, which helps to better focus the reconstruction to an area '
                            ' of interest and further reduce the processing time. Regardless '
                            ' of the width value, the tomogram is always reconstructed around '
                            ' the tilt axis (i.e. the centre of the images in the tilt-series).')

        line = group.addLine('Slices in the Y-axis',
                             expertLevel=LEVEL_ADVANCED,
                             help=' By default, the output tomogram has as many slices along '
                                  ' the tilt axis as the Y dimension of the images in the tilt-series. '
                                  ' This option allows the user to specify the set of slices to be '
                                  ' reconstructed. Thus, this option helps to better focus the '
                                  ' reconstruction to an area of interest and further reduce the '
                                  ' processing time. The user has to specify two values y1,y2 numbered '
                                  ' from 0 (i.e. in the range [0,Ny-1], with Ny being the Y dimension '
                                  ' of the images), which represent the indices of the first and last '
                                  ' slices to be included in the output tomogram.')

        line.addParam('iniSlice', IntParam, default=0, label='Initial')
        line.addParam('finSlice', IntParam, default=0, label='Final')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        pass

    def rotXTomo(self, tsId):
        """Result of the reconstruction must be rotated 90 degrees around the X
        axis to recover the original orientation (due to jjsoft design)"""
        inTomoFile = self._getTmpTomoOutFName(tsId)
        outTomoFile = self._getTomoOutFName(tsId)

        with mrcfile.open(inTomoFile, mode='r', permissive=True) as mrc:
            rotData = np.rot90(mrc.data)

        with mrcfile.open(outTomoFile, mode='w+') as mrc:
            mrc.set_data(rotData)

    def createOutputStep(self, tsId):
        ts = self.tsDict[tsId]
        acq = ts.getAcquisition().clone()
        outputTomos = self._getOutputSetOfTomograms()
        tomo = Tomogram()
        tomo.setTsId(ts.getTsId())
        tomo.setLocation(self._getTomoOutFName(tsId))
        tomo.setSamplingRate(ts.getSamplingRate())
        tomo.setAcquisition(acq)
        tomo.setOrigin()
        outputTomos.append(tomo)
        outputTomos.update(tomo)
        outputTomos.write()
        self._store(outputTomos)

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        errorMsg = []
        if self.height.get() % 2 == 1:
            errorMsg.append('The thickness must be an even number')

        return errorMsg


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

    def getWorkingDirName(self, tsId):
        return self._getTmpPath(tsId)

    def _getOutputSetOfTomograms(self):
        outTomograms =  getattr(self, outputTomoRecObjects.tomograms.name, None)
        if outTomograms:
            outTomograms.enableAppend()
            tomograms = outTomograms
        else:
            inTsSet = self.inputSetOfTiltSeries.get()
            tomograms = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            tomograms.copyInfo(inTsSet)
            tomograms.setStreamState(Set.STREAM_OPEN)
            setattr(self, outputTomoRecObjects.tomograms.name, tomograms)
            self._defineOutputs(**{self._OUTNAME: tomograms})
            self._defineSourceRelation(inTsSet, tomograms)

        return tomograms

    def _getTomoOutFName(self, tsId):
        return join(self._getTsExtraDir(tsId), '%s.mrc' % tsId)

    def _getTmpTomoOutFName(self, tsId):
        return join(self._getTsTmpDir(tsId), '%s.mrc' % tsId)

    def _getTsTmpDir(self, tsId):
        return self._getTmpPath(tsId)

    def _getTsExtraDir(self, tsId):
        return self._getExtraPath(tsId)


