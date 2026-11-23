from Barwy import Spectrum
import numpy as np
import matplotlib.pyplot as plt

def main():
    fig, (ax1,ax2) = plt.subplots(1,2,layout="tight", figsize=(12, 4))

    headerlines = 19
    x = []
    y = []
    with open("LED_zolty.txt", 'r') as file:
        for n in range(headerlines): file.readline()
        for line in file:
            line = line.split('	')
            line = [line[n].strip(' \n') for n in range(len(line))]
            x_new = np.float64(line[0])
            y_new = np.float64(line[1])
            x.append(x_new)
            y.append(y_new)

    LED = Spectrum(x=x, spectrum=y, illuminant='d65', normalize=True)
    LED.plot_functions()
    ax2.fill_between([0, 1], 1, facecolor=LED.get_rgb())
    LED.plot_spectrum(ax=ax1)
    plt.show()

if __name__ == "__main__":
    main()
