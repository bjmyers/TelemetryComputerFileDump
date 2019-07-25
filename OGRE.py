import numpy as np
from prtp.Sources import Subannulus
from prtp.WolterOptic import WolterModule
from prtp.GratingStack import GratingStack
import astropy.units as u


class OgreMirrorModule:
    '''
    Class OgreMirrorModule:
    Stores the parameters for the Mirror Module for OGRE
    
    To get the WolterModule object, use the following code:
    >> from prtp.Missions.OGRE import OgreMirrorModule
    >> mirror = OgreMirrorModule()
    '''
    def __new__(cls):
        r0s = np.array([165., 167.5503, 170.1193, 172.7023,
                        175.3143, 177.9404, 180.5859, 183.2509,
                        185.9355, 188.6398, 191.3640, 194.1083]) * u.mm
        
        z0s = np.ones(12) * 3500. * u.mm
        ax_lens = np.ones(12) * 100 * u.mm
        mir_seps = np.ones(12) * 5 * u.mm
        
        return WolterModule(r0=r0s,z0=z0s,axial_length=ax_lens,mirror_sep=mir_seps,beckmann_scatter=True,ripple=1.48e-5)

#TODO: OgreSource in progress, see Ben for what he wants
class OgreSource:
    '''
    Class OgreSource:
    Stores the parameters for the Source for OGRE
    
    To get the Source object, use the following code:
    >> from prtp.Missions.OGRE import OgreSource
    >> source = OgreSource()
    >> source.num = num_photons
    This last line allows you to change the number of photons that will be emitted
    '''
    def  __new__(cls):
        z0 = 3500.
        mirror_length = 100.
        mirror_sep = 5.
        
        min_radial_extent = conic.primrad(z0 + mirror_sep/2, r_int[0], z0) * u.mm
        max_radial_extent = conic.primrad(z0 + mirror_length + mirror_sep/2,
                                            r_int[-1], z0) * u.mm
        
        return Subannulus(num=1000,min_radial_extent,max_radial_extent,
                            dphi=60*u.deg,z=-z0)
        # wave and order are currently None