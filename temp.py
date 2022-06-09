from numpy import zeros
from cmath import exp, pi
from numpy import loadtxt
from pylab import plot, xlim, show
import numpy as np

# discrete fourier transform 

def dft (x):
    print("in dft")
    N = len (x)
    c = zeros (N//2 + 1, complex)
    for k in range (N//2 + 1):
        for n in range (N):
            c[k] += x[n] * exp (-2j * pi * k * n / N)
    return c 

# recursive solution to FFT
def fft (x):
    x = np.asarray (x, dtype = float)
    N = x.shape[0]
    if (N % 2 > 0):
        raise ValueError ("must be a power of 2")
    elif (N <= 2):
        return dft (x)
    else:
        X_even = fft (x[::2])
        X_odd = fft (x[1::2])
        terms = np.exp (-2j * np.pi * np.arange(N) / N)
        return np.concatenate(X_even + terms[:int(N/2)] * X_odd, X_even + terms[int(N/2):] * X_odd)

def main():
    x = loadtxt ("./trumpet.txt", float)
    print(x, len(x))
    x = x[0:2048]
    #c = np.fft.fft (x)
    #print(c)
    d = fft(x)
    plot (abs(d))
    xlim (0, len(x)/2)
    show()
    return()
main()