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
from .constants import JJSOFT_HOME

__version__ = '3.0.2'


class Plugin(pwem.Plugin):
    _homeVar = JJSOFT_HOME
    _pathVars = [JJSOFT_HOME]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(JJSOFT_HOME, "jjsoft")

    @classmethod
    def getTomoBFlowProgram(cls):
        return join(cls.getHome(), 'tomobflow', 'bin', 'tomobflow')

    @classmethod
    def getTomoEEDProgram(cls):
        return join(cls.getHome(), 'tomoeed', 'bin', 'tomoeed')

    @classmethod
    def getAlignProgramsPath(cls):
        return join(cls.getHome(), 'tomoalign_Jan2019_linux', 'bin')

    @classmethod
    def getTomowarpalignProgram(cls):
        return join(cls.getAlignProgramsPath(), 'tomowarpalign')

    @classmethod
    def getTomoAlignProgram(cls):
        return join(cls.getAlignProgramsPath(), 'tomoalign')

    @classmethod
    def getTomoRecProgram(cls):
        return join(cls.getHome(), 'tomo3d')

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith("")

