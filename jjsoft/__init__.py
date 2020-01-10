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
This package contains the protocols and data for jjsoft.0
"""
import os
import pyworkflow.em

from pyworkflow.utils import Environ
from .constants import JJSOFT_HOME

class Plugin(pyworkflow.em.Plugin):
    _homeVar = JJSOFT_HOME
    _pathVars = [JJSOFT_HOME]


    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(JJSOFT_HOME, "tomobins")


    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Xmipp2. """
        environ = Environ(os.environ)

        environ.update({
            'PATH': os.path.join(cls.getVar(JJSOFT_HOME))
        }, position=Environ.BEGIN)

        return environ

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith("")


pyworkflow.em.Domain.registerPlugin(__name__)


@classmethod
def defineBinaries(cls, env):
    """ Define the Xmipp binaries/source available tgz.
    """
    scons = tryAddPipModule(env, 'scons', '3.0.4', default=True)

def tryAddPipModule(env, moduleName, *args, **kwargs):
    """ To try to add certain pipModule.
        If it fails due to it is already add by other plugin or Scipion,
          just returns its name to use it as a dependency.
        Raise the exception if unknown error is gotten.
    """
    try:
        return env.addPipModule(moduleName, *args, **kwargs)._name
    except Exception as e:
        if "Duplicated target '%s'" % moduleName == str(e):
            return moduleName
        else:
            raise Exception(e)