# **************************************************************************
# *
# * Authors:  Laura del Cano (ldelcano@cnb.csic.es)
# *           Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es) [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
"""
This package contains the protocols and data for tomo3d and related software
"""
from os.path import join
import pwem
from .constants import *
__version__ = '3.1.3'


class Plugin(pwem.Plugin):
    _homeVar = TOMO3D_HOME_VAR
    _pathVars = [TOMO3D_HOME_VAR]
    _url = 'https://github.com/scipion-em/scipion-em-tomo3d'

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(TOMO3D_HOME_VAR, TOMO3D_DEFAULT_HOME)

    @classmethod
    def getTomoBFlowProgram(cls):
        return join(cls.getHome(), 'tomobflow')

    @classmethod
    def getTomoEEDProgram(cls):
        return join(cls.getHome(), 'tomoeed')

    @classmethod
    def getTomowarpalignProgram(cls):
        return join(cls.getHome(), 'tomowarpalign')

    @classmethod
    def getTomoAlignProgram(cls):
        return join(cls.getHome(), 'tomoalign')

    @classmethod
    def getTomo3dProgram(cls):
        return join(cls.getHome(), TOMO3D)

    @classmethod
    def getTomoRecProgram(cls):
        return join(cls.getHome(), 'tomorec')

    @classmethod
    def defineBinaries(cls, env):
        TOMO3D_INSTALLED = '%s.txt' % TOMO3D
        dlZipFile = TOMO3D + '.zip'
        installationCmd = 'wget %s -O %s && ' % (TOMO3D_BIN_URL, dlZipFile)
        installationCmd += 'unzip %s -d %s && ' % (dlZipFile, cls.getHome())
        installationCmd += 'touch %s' % TOMO3D_INSTALLED  # Flag installation finished

        env.addPackage(TOMO3D,
                       version=TOMO3D_VERSION,
                       tar='void.tgz',
                       commands=[(installationCmd, TOMO3D_INSTALLED)],
                       neededProgs=["wget", "tar", "zip"],
                       default=True)


