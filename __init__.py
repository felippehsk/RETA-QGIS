"""
RETA - Robust Error-treatment Toolkit for Agriculture

QGIS plugin for spatial data filtering of agricultural sensor data.
Adapted from Yield Editor 2.0 (Sudduth, Drummond & Myers, 2012).
"""


def classFactory(iface):
    from .reta_plugin import RetaPlugin
    return RetaPlugin(iface)
