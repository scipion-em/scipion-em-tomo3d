# -*- coding: utf-8 -*-
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
"""
This module implement some wizards
"""
from pyworkflow.wizard import Wizard
from .protocols import ProtJjsoftProtDenoiseTomogram
from .protocols.protocol_denoise_tomogram import DENOISE_EED


class Tomo3dDenoiseIterationsWizard(Wizard):
    _targets = [(ProtJjsoftProtDenoiseTomogram, ['nIter'])]

    def _getNiter(self, protocol):
        if protocol.method.get() == DENOISE_EED:
            nIter = 10
        else:  # BFlow
            nIter = 70
        return nIter

    def show(self, form):
        form.setVar('nIter', self._getNiter(form.protocol))
