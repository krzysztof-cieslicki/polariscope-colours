import matplotlib.pyplot as plt
import numpy as np
from colormath.color_conversions import convert_color
from colormath.color_objects import sRGBColor, xyYColor

class Spectrum:

    def __init__(self, filename=None, x=None, spectrum=None, illuminant='c', headerlines=None, delimiter="\t"):
        if filename is not None:
            x = []
            spectrum = []
            peak = 0
            with open(filename, 'r') as file:
                for n in range(headerlines): file.readline()
                for line in file:
                    line = line.split(delimiter)
                    line = [float(line[n].strip(' \n')) for n in range(len(line))]
                    if line[1] > peak: peak = line[1]
                    x.append(line[0])
                    spectrum.append(line[1])

        self.__lambda = x
        self.__spec = np.array(spectrum)/max(spectrum)*1024
        self.__illuminant = illuminant

        self.__x_bar = self.__g(x,(1.056, 599.8, 0.0264, 0.0323), (0.362, 442, 0.0624, 0.0374),
                                (-0.065, 501.1, 0.0490, 0.0382))
        self.__y_bar = self.__g(x, (0.821, 568.8, 0.0231, 0.0247), (0.286, 530.9, 0.0613, 0.0322))
        self.__z_bar = self.__g(x, (1.217, 437.0, 0.0845, 0.0278), (0.681, 459.0, 0.0385, 0.0725))

        x_scale = 100/np.trapezoid(self.__x_bar, x=self.__lambda)
        y_scale = 100/np.trapezoid(self.__y_bar, x=self.__lambda)
        z_scale = 100/np.trapezoid(self.__z_bar, x=self.__lambda)

        self.__x_bar *= x_scale
        self.__y_bar *= y_scale
        self.__z_bar *= z_scale

    def __g(self, x, *coeffs):
        val = np.zeros(len(x))
        for line in coeffs:
            A, mu, t1, t2 = line
            for n in range(len(x)):
                arg = x[n]
                if arg < mu:
                    new_val = np.exp(-(t1 * (arg - mu)) ** 2 / 2)
                else:
                    new_val = np.exp(-(t2 * (arg - mu)) ** 2 / 2)

                new_val *= A
                val[n] += new_val

        return val

    def plot_spectrum(self, ax=None, xlim=(380, 800)):
        if ax is None:
            fig, ax = plt.subplots(layout='tight', figsize=(8,6))

        y_max = min(self.__spec)
        for n in range(len(self.__lambda)):
            arg = self.__lambda[n]
            if arg >= xlim[0]:
                if arg <= xlim[1]:
                    val = self.__spec[n]
                    if val > y_max:
                        y_max = val
        y_max = np.round(y_max*10, 0)/10
        ax.plot(self.__lambda, self.__spec, color='r', linewidth=1)
        ax.set_xlabel('Wavelength [nm]', fontsize=16)
        ax.set_ylabel('Intensity [a.u.]', fontsize=16)

        ax.grid(which='major', visible=True, linewidth=1)
        ax.grid(which='minor', visible=True, linewidth=0.5, linestyle='dotted')
        ax.minorticks_on()
        ax.tick_params(axis='both', labelsize=16)
        ax.set_xlim(xlim)
        ax.set_ylim((0, y_max*1.01))

        return ax

    def plot_functions(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots(layout='tight')

        ax.plot(self.__lambda, self.__x_bar, color='r')
        ax.plot(self.__lambda, self.__y_bar, color='g')
        ax.plot(self.__lambda, self.__z_bar, color='b')
        ax.set_xlabel(r'$\lambda\:[\rm nm]$')
        ax.legend([r'$\overline{x}(\lambda)$', r'$\overline{y}(\lambda)$', r'$\overline{z}(\lambda)$'], fancybox=False,
                   fontsize=13)
        ax.set_xlim([380, 780])
        ax.set_ylim([0, 2])
        print(f'x:{np.trapezoid(self.__x_bar, x=self.__lambda)}\ny:{np.trapezoid(self.__y_bar, x=self.__lambda)}\nz:{np.trapezoid(self.__z_bar, x=self.__lambda)}')
        return ax

    def get_xycolor(self):

        X = np.trapezoid(self.__spec * self.__x_bar, x=self.__lambda)
        Y = np.trapezoid(self.__spec * self.__y_bar, x=self.__lambda)
        Z = np.trapezoid(self.__spec * self.__z_bar, x=self.__lambda)

        #X, Y, Z = X * self.__scale[0], Y * self.__scale[1], Z * self.__scale[2]
        pass
        if X+Y+Z < 1e-12: x, y, z = (0,0,0)
        else:
            x = X / (X+Y+Z)
            y = Y / (X+Y+Z)
            z = Z / (X+Y+Z)
        return x, y, Y

    def get_rgb(self):
        x, y, Y = self.get_xycolor()
        xyY = xyYColor(x, y, Y, illuminant=self.__illuminant)
        RGB_converted = convert_color(xyY, sRGBColor, target_illuminant=self.__illuminant)
        r,g,b = RGB_converted.get_value_tuple()
        RGB = sRGBColor(r,g,b, is_upscaled=True)
        return RGB.get_value_tuple()

    def updated_spec(self, transmission):
        return self.__spec * transmission

    def get_functions(self):
        return self.__x_bar, self.__y_bar, self.__z_bar

    def get_lambda(self):
        return self.__lambda

    def get_spec(self):
        return self.__spec

    def get_illum(self):
        return self.__illuminant

class Polariscope:

    def __init__(self, crossed=True):
        self.__crossed = crossed

    def get_transmission(self, x, opd=None):
        gamma = np.pi*opd/x
        transmitted = np.sin(gamma) ** 2
        if self.__crossed is False:
            transmitted = 1 - transmitted
        return np.array(transmitted)

class Crystal:

    def __init__(self, d, n1, n2, *, name=None):
        self.__d = d # ile w mm
        self.__n1 = n1
        self.__n2 = n2
        self.__name = name

    def get_biref(self, x=555, nm=True):
        if nm is True:
            x = np.multiply(x, 1e-3)

        return self.__n2(x)-self.__n1(x)

    def get_opd(self, wavelen, nm=True):
        opd = self.get_biref(x=wavelen, nm=nm) * self.__d / wavelen
        return opd

def draw_outline(ax=None):
    file = open("CIE_cc_1931_2deg.csv", 'r')
    x = []
    y = []
    for line in file:
        line = line.split(",")
        line = [np.float64(line[n].strip(' \n')) for n in range(0, len(line))]
        _ , new_x, new_y, _ = line
        x.append(new_x)
        y.append(new_y)

    file.close()
    if ax is None:
        fig, ax = plt.subplots(layout='tight')

    ax.plot(x, y, linewidth=0.75, color='k')
    ax.plot([x[0], x[len(x)-1]], [y[0], y[len(y)-1]], linewidth=0.75, color='k')
    ax.grid(which='major', visible=True, linewidth=1)
    ax.grid(which='minor', visible=True, linewidth=0.5, linestyle='dotted')
    ax.minorticks_on()
    ax.tick_params(axis='both', labelsize=16)
    ax.set_xlabel('x', fontsize=16)
    ax.set_ylabel('y', fontsize=16)
    ax.set_aspect('equal')
    return ax

def draw_srgb(ax=None):
    if ax is None:
        fig, ax = plt.subplot()

    ax.plot([0.64, 0.30, 0.15, 0.64], [0.33, 0.60, 0.06, 0.33], linewidth=0.5, color='k')

def color_curve(x, light, opd_min=0, opd_max=2000, opd_points=2000, crossed=True):
    fig1, ax1 = plt.subplots(figsize=(8, 6), layout='tight')

    polaryskop = Polariscope(crossed)
    if crossed is True:
        if opd_min == 0:
            opd_min = 1e-3

    opd = np.linspace(opd_min, opd_max, opd_points)
    x_plot = []
    y_plot = []
    rgb = []
    peak = 1
    for n in range(opd_points):
        new_spectrum = light.get_spec() * polaryskop.get_transmission(x, opd[n])
        exiting = Spectrum(x=x, spectrum=new_spectrum, illuminant=light.get_illum())
        xy_new = exiting.get_xycolor()
        x_plot.append(xy_new[0])
        y_plot.append(xy_new[1])
        r, g, b = exiting.get_rgb()
        #if max(r,g,b) > peak: peak = max(r,g,b)
        rgb.append((r, g, b))

    rgb = np.array(rgb)/peak

    draw_outline(ax1)
    draw_srgb(ax1)
    ax1.plot(x_plot, y_plot, color='r', linewidth=1)
    dx = opd[1] - opd[0]

    fig2, ax2 = plt.subplots(layout='tight')
    for n in range(opd_points):
        ax2.fill_between((opd[n], opd[n] + dx), opd_max * 0.1, color=rgb[n])
    ax2.set_aspect('equal')
    ax2.set_xlabel(r'$\Delta L\:[\rm nm]$', fontsize=16)
    ax2.get_yaxis().set_visible(False)
    ax2.set_xlim([0, opd_max])
    ax2.set_ylim([0, opd_max*0.1])
    ax2.tick_params(labelsize=16)


def main():

    illuminant = 'CIE_illum_C.csv'

    incident = Spectrum(filename=f"./{illuminant}", illuminant='C', headerlines=1, delimiter=",")
    x = incident.get_lambda()
    incident.plot_functions()
    #fig1.savefig('color matching functions', dpi=300)

    color_curve(x=x, opd_min=0, opd_max=2000, opd_points=1000, crossed=True, light=incident)
    plt.show()

if __name__ == "__main__":
    main()