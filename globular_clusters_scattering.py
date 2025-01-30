"""A python file to use pygedm to calculate the scattering
timescales for clusters

@author: Deven Bhakta
@date: sent on 1/28/25
@title: globular_clusters_scattering.py
@CPU: MacBook Pro Apple M3
@Operating System: Sonoma 14.6.1
@Interpreter and version no.: python 3.12.2
"""

# %%
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

from astroplan.plots import plot_altitude, plot_sky
from astroplan import FixedTarget, Observer, AltitudeConstraint, Target
from astropy.time import Time
from astropy.coordinates import SkyCoord

from mwprop.scattering_functions.scattering_functions2020 import tauiss
from mwprop.ne2001p.NE2001 import ne2001

from pygedm import dist_to_dm

# from presto import psr_utils as psu

# %%
source_list = ["Liller 1",'NGC 6388', 'NGC 6139',"NGC 6304","NGC 6717","Terzan 2"] # 2MS-CO1 and GLIMPSE-CO2 need to be added manually
targets = []
for source in source_list:
    targets.append(FixedTarget.from_name(source))

coords = SkyCoord(ra='18h08m21.8s', dec = "-19d49m47s",frame='icrs')
targets.append(FixedTarget(coord=coords,name='2MS-CO1'))

coords = SkyCoord(ra='18h18m30.5s', dec = "-16d58m38s",frame='icrs')
targets.append(FixedTarget(coord=coords,name='GLIMPSE-CO2'))
# %%
dist = [8.06,11.15,10.04,5.9,7.1,7.5,3.6,5.5]
freq = 2.165
tau_ne = []
DM_ne = []
tau_ywm = []
DM_ywm = []
for d,target in zip(dist,targets):
    Dk,Dv,Du,Dd=ne2001(ldeg=target.coord.galactic.l.degree,bdeg=target.coord.galactic.b.degree,dmd=d,ndir=-1,dmd_only=False)
    tau_ne.append(tauiss(d,Dv['SM'],freq)*u.ms)
    DM_ne.append(Dv['DM']*u.pc/(u.cm**3))
    tau_ywm.append(dist_to_dm(target.coord.galactic.l.degree, target.coord.galactic.b.degree, dist = d*1000,nu=freq)[1])
    DM_ywm.append(dist_to_dm(target.coord.galactic.l.degree, target.coord.galactic.b.degree, dist = d*1000,nu=freq)[0])
# %%
print(f'Freq: {freq} GHz')
for i,target in enumerate(targets):
    print(f'{target.name} \t NE2001: tau - {tau_ne[i].to(u.us):.1f}| DM - {DM_ne[i]:.2f} \t YMW16: tau - {tau_ywm[i].to(u.us):.1f}| DM - {DM_ywm[i]:.2f}')
    # print(f'The NE2001 scattering timescale for {target.name} at {freq} GHz: \t {tau_ne[i].to(u.us):.3f} \t DM: {DM_ne[i]:.2f}' )
    # print(f'The YWN16 scattering timescale for {target.name} at {freq} GHz: \t {tau_ywm[i].to(u.us):.3f} \t DM: {DM_ywm[i]:.2f}')
# %%
dist_to_dm(targets[0].coord.galactic.l.degree, targets[0].coord.galactic.b.degree, dist = dist[0]*1000,nu=freq)[1]

# %%
for target in targets:
    print(f'{target.name} \t {target.coord.galactic.l.degree} \t {target.coord.galactic.b.degree}')
# %%
