"""
This program will simulate surface plasmon resonance of MoS2 phototransistor. 
We shall create the geometry where MoS2 is placed over SiO2 and Si substrate.
Silver nanodisks (AgND) will be placed on top of MoS2 layer.

Created on Sun Jan  7 10:48:38 2024

@author: raonaq
"""
import math
import matplotlib.pyplot as plt
import numpy as np
import meep as mp
import time

# Import materials
from meep.materials import Si
from meep.materials import Ag
from meep.materials import SiO2

start_time = time.time()

resolution = 60  # pixels/μm

# dpml = 1.0  # PML thickness
# dsub = 3.0  # substrate thickness
# dpad = 3.0  # padding between grating and PML
# gp = 10.0  # grating period
# gh = 0.5  # grating height
# gdc = 0.5  # grating duty cycle

w = 1 # Overall width of the device
dSi = 0.6 # Thickness of Si substrate
dSiO2 = 0.3 # Thickness of SiO2 layer
dMoS2 = 0.0007 # Thickness of MoS2 layer
dAg = 0.04 # Thickness of AgND
wAg = 0.160 # Width of AgND
pAg = 0.260 # Period of AgND
dpad = 0.2 # Padding between AgND and PML
dpml = 0.2 # PML thickness 

sy = dpml+dSi+dSiO2+dMoS2+dAg+dpad+dpml # Cell size along y-axis
sx = pAg # Cell size along x-axis

cell_size = mp.Vector3(sx, sy, 0)
pml_layers = [mp.PML(thickness=dpml, direction=mp.Y)]

wvl_min = 0.4  # min wavelength
wvl_max = 0.6  # max wavelength
fmin = 1 / wvl_max  # min frequency
fmax = 1 / wvl_min  # max frequency
fcen = 0.5 * (fmin + fmax)  # center frequency
df = fmax - fmin  # frequency width

# Add Source
# src_pt = mp.Vector3(0,-0.5 * sx + dpml + 0.5 * dsub,0) # Source on bottom
src_pt = mp.Vector3(0,0.5*sy-dpml-0.5*dpad,0) # Source on top

# Pulsed source
# sources = [
#     mp.Source(
#         mp.GaussianSource(fcen, fwidth=df),
#         component=mp.Ez,
#         center=src_pt,
#         size=mp.Vector3(sx,0,0),
#     )
# ]

# Continuous source
sources = [
    mp.Source(
        mp.ContinuousSource(frequency=1.89), # wl=0.532 um converted to f=1/wl 
        component=mp.Ez,
        center=src_pt,
        size=mp.Vector3(sx,0,0)
    )
]

k_point = mp.Vector3(0, 0, 0)

glass = mp.Medium(index=1.5)

symmetries = [mp.Mirror(mp.X)]

sim = mp.Simulation(
    resolution=resolution,
    cell_size=cell_size,
    boundary_layers=pml_layers,
    k_point=k_point,
    # default_material=glass,
    sources=sources,
    symmetries=symmetries,
)

nfreq = 21

# Add monitor
# mon_pt = mp.Vector3(0,0.5 * sx - dpml - 0.5 * dpad,0) # Monitor on top
# mon_pt = mp.Vector3(0,-0.5 * sx + dpml + 0.5 * dsub,0) # Monitor on bottom
mon_pt = mp.Vector3()  # One point monitor in the middle

flux_mon = sim.add_flux(
    fcen, df, nfreq, mp.FluxRegion(center=mon_pt, size=mp.Vector3(sx,0,0))
)

f = plt.figure(dpi=120)
sim.plot2D(ax=f.gca())
plt.show()

# sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mon_pt, 1e-9)) # Until decayed
sim.run(until=200) # Until a specified time step

input_flux = mp.get_fluxes(flux_mon)

sim.reset_meep()

geometry = [
    
    # Si substrate
    mp.Block(
        material=Si,
        size=mp.Vector3(w,dpml+dSi,w),
        center=mp.Vector3(0,-0.5*sy+0.5*(dpml+dSi),0),
    ),
    
    # SiO2 layer
    mp.Block(
        material=SiO2,
        size=mp.Vector3(w,dSiO2,w),
        center=mp.Vector3(0,-0.5*sy+(dpml+dSi)+0.5*dSiO2,0),
    ),
    
    # MoS2 layer
    mp.Block(
        material=mp.Medium(epsilon=12.46), # n = 3.53 at 530 um
        size=mp.Vector3(w,dMoS2,w),
        center=mp.Vector3(0,-0.5*sy+(dpml+dSi+dSiO2)+0.5*dMoS2,0),
    ),
    
    # AgND
    mp.Block(
        material=Ag,
        size=mp.Vector3(wAg,dAg,wAg),
        center=mp.Vector3(0,-0.5*sy+(dpml+dSi+dSiO2+dMoS2)+0.5*dAg,0),
    ),
    
]

sim = mp.Simulation(
    resolution=resolution,
    cell_size=cell_size,
    boundary_layers=pml_layers,
    geometry=geometry,
    k_point=k_point,
    sources=sources,
    symmetries=symmetries,
)

# mode_mon = sim.add_flux(
#     fcen, df, nfreq, mp.FluxRegion(center=mon_pt, size=mp.Vector3(gdc*gp,0,0))
# )

f2 = plt.figure(dpi=600)
sim.plot2D(ax=f2.gca())
plt.show()

# sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mon_pt, 1e-9)) # Until decayed
sim.run(until=200) # Until 200 time steps

# Plot fields
plt.figure(dpi=600)
sim.plot2D(fields=mp.Ez)
plt.show()

end_time = time.time()
elapsed_time = end_time - start_time


