from crystals import quartz, sapphire
from Barwy import Crystal
import matplotlib.pyplot as plt
from numpy import linspace

fig, ax = plt.subplots(layout="tight", figsize=(6,4))
sapphire = Crystal(n1=sapphire["n1"], n2=sapphire["n2"], name="Sapphire", d=None)
quartz = Crystal(n1=quartz["n1"], n2=quartz["n2"], name="Quartz", d=None)

waves = linspace(.4, .8)
ax.plot(waves, sapphire.get_biref(waves), color="C3", linewidth=2, label="Sapphire")
ax.plot(waves, quartz.get_biref(waves), color="C2", linewidth=2, label="Quartz")
ax.set_xlabel(r"Wavelength ($\mu$m)", fontsize=16)
ax.set_ylabel(r"Birefringence", fontsize=16)
ax.ticklabel_format(axis="y", scilimits=(-2, 2))
ax.tick_params(axis="both", labelsize=14)
ax.legend(fontsize=14)
plt.show()


