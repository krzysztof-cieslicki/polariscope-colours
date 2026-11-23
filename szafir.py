from Barwy import Crystal, Spectrum
import numpy as np
import matplotlib.pyplot as plt


thickness = 38 # [um]
fig, ax = plt.subplots(layout='tight', figsize=(7.5,5))
fig2, ax2 = plt.subplots(layout='tight', figsize=(4,4))

sapphire = Crystal(n1=lambda x: (1+1.4313493/(1-(0.0726631/x)**2)+0.65054713/(1-(0.1193242/x)**2)+5.3414021/(1-(18.028251/x)**2))**.5,
                   n2=lambda x: (1+1.5039759/(1-(0.0740288/x)**2)+0.55069141/(1-(0.1216529/x)**2)+6.5927379/(1-(20.072248/x)**2))**.5,
                   d=thickness * 1e3)

light = Spectrum(illuminant='D65', filename="illuminant files/CIE_std_illum_D65.csv", headerlines=0, delimiter=',', normalize=True)
wavelength = light.get_lambda()

transmission = np.sin(np.pi * sapphire.get_opd(wavelength*1e-3) / wavelength) ** 2
modified = Spectrum(x=wavelength, spectrum=light.updated_spec(transmission))

ax.plot(wavelength, light.get_spec(), linestyle='--', alpha=0.75, color='k', label='Original spectrum')
ax.plot(wavelength, modified.get_spec(), color='r', label='Transmitted spectrum')
ax.set_xlabel("Wavelength [nm]", fontweight="bold", fontsize=16)
ax.set_ylabel("Intensity [a.u.]", fontweight="bold", fontsize=16)
ax.set_title(f"Sapphire sample\n{thickness} um thickness", fontweight="bold", fontsize=18)
ax.legend(fontsize=14)
ax.tick_params(axis="both", labelsize=14)
ax.ticklabel_format(axis="y", scilimits=(-2, 2))

ax2.fill_between(x=(0, 1), y1=1, y2=0, color=modified.get_rgb())
ax2.xaxis.set_visible(False)
ax2.yaxis.set_visible(False)
plt.show()
