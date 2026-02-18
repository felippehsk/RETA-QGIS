"""
RETA Processing Provider for QGIS.

Registers the Spatial Data Filtering algorithm in the QGIS Processing Toolbox
under the RETA provider group.
"""

from qgis.core import QgsProcessingProvider
from .spatial_data_filtering import SpatialDataFilteringAlgorithm


class RetaProvider(QgsProcessingProvider):

    def loadAlgorithms(self):
        self.addAlgorithm(SpatialDataFilteringAlgorithm())

    def id(self):
        return 'reta'

    def name(self):
        return 'RETA'

    def longName(self):
        return 'RETA - Robust Error-treatment Toolkit for Agriculture'
