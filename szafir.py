from Barwy import Crystal, Spectrum
import numpy as np
import matplotlib.pyplot as plt
from crystals import sapphire

thickness = 80 # [um]
fig, ax = plt.subplots(layout='tight', figsize=(7.5,5))
fig2, ax2 = plt.subplots(layout='tight', figsize=(4,4))

sapphire = Crystal(n1=sapphire["n1"], n2=sapphire["n2"], d=thickness * 1e3)
light = Spectrum(illuminant='D65', filename="illuminant files/CIE_std_illum_A_1nm.csv", headerlines=0, delimiter=',', normalize=True)
wavelength = light.get_lambda()

transmission = np.sin(np.pi * sapphire.get_opd(wavelength*1e-3) / wavelength) ** 2
modified = Spectrum(x=wavelength, spectrum=light.get_spec()*transmission)

ax.plot(wavelength, light.get_spec(), linestyle='--', alpha=0.75, color='k', label='Original spectrum')
ax.plot(wavelength, modified.get_spec(), color='r', label='Transmitted spectrum')
ax.set_xlabel("Wavelength [nm]", fontsize=16)
ax.set_ylabel("Intensity [a.u.]", fontsize=16)
ax.set_title(f"Sapphire sample\n{thickness} um thickness", fontweight="bold", fontsize=18)
ax.legend(fontsize=14)
ax.tick_params(axis="both", labelsize=14)
ax.ticklabel_format(axis="y", scilimits=(-2, 2))

ax2.fill_between(x=(0, 1), y1=1, y2=0, color=light.modified_spec_rgb(transmission))
ax2.set_title("Simulated colour", fontweight="bold", fontsize=16)
for axis in (ax2.xaxis, ax2.yaxis):
    axis.set_visible(False)

plt.show()
