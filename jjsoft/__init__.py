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
This package contains the protocols and data for jjsoft
"""
from os.path import join
import pwem
from pyworkflow.utils import yellowStr
from .constants import JJSOFT_HOME, JJSOFT, JJSOFT_README, JJSOFT_BIN_VERSION, TOMO3D

__version__ = '3.0.3'


class Plugin(pwem.Plugin):
    _homeVar = JJSOFT_HOME
    _pathVars = [JJSOFT_HOME]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(JJSOFT_HOME, JJSOFT + '-%s' % JJSOFT_BIN_VERSION)

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
        # At this point of the installation execution cls.getHome() is None, so the em path should be provided
        # Only the directory will be generated, because the binaries must be downloaded manually from José
        # Jesús website, filling a form
        JJSOFT_INSTALLED = 'readme.txt'
        msg = 'Binaries must be installed manually. Please follow the instructions described here: ' + JJSOFT_README
        # installationCmd = 'touch %s &&' % JJSOFT_INSTALLED  # Flag installation finished
        installationCmd = 'echo %s && ' % yellowStr(msg.upper())
        installationCmd += 'echo "%s" >> %s' % (msg, JJSOFT_INSTALLED)

        env.addPackage(JJSOFT,
                       version=JJSOFT_BIN_VERSION,
                       tar='void.tgz',
                       commands=[(installationCmd, JJSOFT_INSTALLED)],
                       neededProgs=["wget", "tar"],
                       default=True)

#TODO validateInstallation


