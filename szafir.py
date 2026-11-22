from Barwy import Crystal, Spectrum
import numpy as np
import matplotlib.pyplot as plt

sapphire = Crystal(n1=lambda x: (1+1.4313493/(1-(0.0726631/x)**2)+0.65054713/(1-(0.1193242/x)**2)+5.3414021/(1-(18.028251/x)**2))**.5,
                   n2=lambda x: (1+1.5039759/(1-(0.0740288/x)**2)+0.55069141/(1-(0.1216529/x)**2)+6.5927379/(1-(20.072248/x)**2))**.5,
                   d=50e6) # 5 mm w nm

light = Spectrum(illuminant='D65', filename="CIE_std_illum_D65.csv", headerlines=0, delimiter=',')
#light.plot_functions()

wavelength = light.get_lambda()

transmission = np.sin(np.pi * sapphire.get_opd(wavelength) / wavelength) ** 2
modified = Spectrum(x=wavelength, spectrum=light.updated_spec(transmission))

fig, ax = plt.subplots(layout='tight', figsize=(8,6))
fig2, ax2 = plt.subplots(layout='tight', figsize=(8,6))
ax.plot(wavelength, light.get_spec(), linestyle='--', alpha=0.75, color='k', label='Original spectrum')
ax.plot(wavelength, modified.get_spec(), color='r', label='Transmitted spectrum')
ax2.plot(wavelength, transmission, color="g", label="Transmission")
ax.set_xlabel("Wavelength [nm]")
ax.set_ylabel("Intensity [a.u.]")
ax.legend(fontsize=13)
ax.set_title("Sapphire sample\n5 mm thickness")

fig3, ax3 = plt.subplots(layout='tight', figsize=(4,4))
rgba = tuple([v/255 for v in modified.get_rgb()])
print(rgba)
ax3.fill_between(x=(0, 1), y1=1, y2=0, color=rgba)
plt.show()
