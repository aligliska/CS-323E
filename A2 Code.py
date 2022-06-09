# Alice Gee, Mohammad Aga, Andrew Yang
# eid: ag67642, mba929, ay6764

import matplotlib.pyplot as plt
import numpy as np
import random as rd
import math

#  x = np.linspace(0,4,100)
#  y = x**2 - np.cos(x) - 3*np.sin(x) - 9
#  plt.plot(x, y, linewidth = 2.0)
#  plt.show()

# given function
f = lambda x: x**2 - np.cos(x) - 3*np.sin(x) - 9

def bisection(f, a, b, tol):
    """An implementation of the Bisection method for finding roots of a function.
    
    Arguments:
        f: the lambda function for the function we want to find the roots for
        a: lower bound
        b: upper bound
        tol: the tolerance
        
    Returns:
        An estimate for the root of the function.
    """
    mid = (a + b) / 2
    new_y = f(mid)
    while (abs(new_y) > tol):
        mid = (a + b) / 2
        new_y = f(mid)
        if (new_y > 0):
            b = mid
        elif (new_y < 0):
            a = mid
        else:
            return mid
    return(mid)

tol = 1e-6
a = 0
b = 5

print("bisection", bisection(f, a, b, tol))


def Regula_Falsi (f, a, b, tol):
    # f Python function or lambda function
    # a lower bound
    # b upper bound
    # tol is the tolerance
    m = (f(b) - f(a)) / (b - a)
    c = (-f(a)/m) + a
    count = 0
    while (abs(f(c)) > tol):
        if (f(b)*f(c) < 0):
            a = c
        elif (f(a)*f(c) < 0):
            b = c
        elif f(c) == 0:
            return (f(c), count)
        count += 1 
        m = (f(b) - f(a)) / (b - a)
        c = (-f(a)/m) + a
 
    return (c, count)

method2 = Regula_Falsi (f, a = 0, b = 4, tol = 1*10**(-6))
print(f"The regula falsi method took {method2[1]} times to find the root: {method2[0]}.")

def newton (f, f_prime, x, tol):
    # f is a differentiable function
    # f_prime is the derivative of f (done by hand)
    # x is initial approximation of the root
    # tol is the tolerance
    
    count = 0
    while (abs(f(x)) > tol):
        x = x - (f(x)/f_prime(x)) 
        count += 1
    return (float(x), count)

f_prime = lambda x: 2*x + np.sin(x) - 3*np.cos(x)
x = 2
method3 = newton(f, f_prime, x, tol = 1*10**(-6))
print(f"The newton method took {method3[1]} times to find the root: {method3[0]}.")


def secant (f, a, b, tol):
    # f is differentiable function
    # a and b are close approximations to the root
    # tol is the tolerance
    f_prime = (f(b)-f(a)) / (b-a)
    c = b - (f(b)/f_prime)
    count = 0
    while abs(f(c)) > tol:
        f_prime = (f(b)-f(a)) / (b-a)
        c = b - (f(b)/f_prime)
        a = b
        b = c
        count += 1
    return(c, count)
    
a = rd.uniform(0,4)
b = rd.uniform(0,4)
method4 = secant (f, a, b, tol = 1*10**(-6))
print(f"The secant method took {method4[1]} times to find the root: {method4[0]}.")


x_0 = 1
f = lambda x: x**2 - np.cos(x) - 3*np.sin(x) - 9
f_prime = lambda x: 2*x + np.sin(x) - 3*np.cos(x)
f_double_prime = lambda x: 3*np.sin(x) + np.cos(x) + 2

def taylor_series(x_0, f, f_prime, f_double_prime, tol = 1*10**(-6)):
    x = (-f_prime(x_0) + math.sqrt((f_prime(x_0))**2 -2*(f(x_0))*(f_double_prime(x_0))) + f_double_prime(x_0)* x_0) / f_double_prime(x_0)
    while abs(f(x)) > tol:
        x = (-f_prime(x_0) + math.sqrt((f_prime(x_0))**2 -2*(f(x_0))*(f_double_prime(x_0))) + f_double_prime(x_0)* x_0) / f_double_prime(x_0)
        x_0 = x
    return(x)

print(f"Test equation #1, where f(x) = x^2 - cos(x) - 3sin(x) - 9:", taylor_series(x_0, f, f_prime, f_double_prime))

x_0 = 3
f = lambda x: np.sin(x) + x**2 - 2 * np.log(x) - 5
f_prime = lambda x: 2*x + np.cos(x) - 2/x
f_double_prime = lambda x: -np.sin(x) + 2 + 2/x**2

def taylor_series(x_0, f, f_prime, f_double_prime, tol = 1*10**(-6)):
    x = (-f_prime(x_0) + math.sqrt((f_prime(x_0))**2 -2*(f(x_0))*(f_double_prime(x_0))) + f_double_prime(x_0)* x_0) / f_double_prime(x_0)
    while abs(f(x)) > tol:
        x = (-f_prime(x_0) + math.sqrt((f_prime(x_0))**2 -2*(f(x_0))*(f_double_prime(x_0))) + f_double_prime(x_0)* x_0) / f_double_prime(x_0)
        x_0 = x
    return(x)

print(f"Test equation #2, where f(x) = sin(x) + x^2 - 2 * log(x) - 5:", taylor_series(x_0, f, f_prime, f_double_prime))

x_0 = 3
f = lambda x: x**2 - 2
f_prime = lambda x: 2*x
f_double_prime = lambda x: 2

def taylor_series(x_0, f, f_prime, f_double_prime, tol = 1*10**(-6)):
    x = (-f_prime(x_0) + math.sqrt((f_prime(x_0))**2 -2*(f(x_0))*(f_double_prime(x_0))) + f_double_prime(x_0)* x_0) / f_double_prime(x_0)
    while abs(f(x)) > tol:
        x = (-f_prime(x_0) + math.sqrt((f_prime(x_0))**2 -2*(f(x_0))*(f_double_prime(x_0))) + f_double_prime(x_0)* x_0) / f_double_prime(x_0)
        x_0 = x
    return(x)

print(f"Test equation #3, where f(x) = x^2 - 2:", taylor_series(x_0, f, f_prime, f_double_prime))


m = 1
g = 9.8
k = 0.1
s_0 = 100
f = lambda t: s_0 - (m*g)/k * t + (m**2)*g/(k**2)(1 - np.exp((-k*t)/m))


