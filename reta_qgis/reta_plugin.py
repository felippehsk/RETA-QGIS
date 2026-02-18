"""
RETA plugin class: registers the Processing provider with QGIS.
"""

from qgis.core import QgsApplication
from .reta_provider import RetaProvider


class RetaPlugin:

    def __init__(self, iface):
        self.iface = iface
        self.provider = None

    def initProcessing(self):
        self.provider = RetaProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):
        self.initProcessing()

    def unload(self):
        if self.provider:
            QgsApplication.processingRegistry().removeProvider(self.provider)
