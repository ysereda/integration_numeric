#!/usr/bin/env python
# coding: utf-8

# # Gauss-Legendre quadrature

# https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature

# In[122]:


import math
def Gauss_Legendre_quadrature(f,a,b,o): # integrand f(x), lower and upper integration limits
    w=np.zeros(o+1); x=np.zeros(o+1); X=np.zeros(o+1);
    if o==1: # 1-point Gaussian
        return (b-a)*f((a+b)/2)
    elif o==2: # 2-point Gaussian
        w[1]=1; w[2]=1; x[2]=1/math.sqrt(3); x[1]= -x[2];
        I2=0;
        for i in range(1,o+1):
            X[i]=(x[i]+1+a)*(b-a)/2;
            I2=I2+w[i]*f(X[i]);
        I2=(b-a)/2*I2;
        return I2;
    elif o==3:
        w[1]= 5/9; w[2]= 8/9; w[3]= 5/9; x[1]=-math.sqrt(3/5); x[2]=0; x[3]=math.sqrt(3/5);
        I3=0;
        for i in range(1, o+1):
            X[i]=(x[i]+1+a)*(b-a)/2;
            I3=I3+w[i]*f(X[i]);
        I3=(b-a)/2*I3;
        return I3;
    elif o==4:
        w[4]= 1/2-1/36*30**(1/2); w[3]= 1/2+1/36*30**(1/2); w[2]=w[3]; w[1]=w[4];
        x[3]= math.sqrt(3/7-2/7*math.sqrt(6/5)); x[4]= math.sqrt(3/7+2/7*math.sqrt(6/5)); x[2]=-x[3]; x[1]=-x[4];
        I4=0;
        for i in range(1, o+1):
            X[i]=(x[i]+1+a)*(b-a)/2;
            I4=I4+w[i]*f(X[i]);
        I4=(b-a)/2*I4;
        return I4;


# In[80]:


# Test function
a=0; b=2.4; # integration interval
def f(x):
    return 2*x/(x**2+1); # integrand


# In[123]:


n_max=4; # maximum order of approximation
I=[0]*(n_max+1); E=[0]*(n_max+1);
I[0]=1.911022890054872722905456; # Exact value of I = int(f(x),x=a..b)
for n in range(1, n_max+1):
    I[n] = Gauss_Legendre_quadrature(f,a,b,n);
    E[n] = abs(I[n]/I[0]-1);


# In[127]:


import numpy as np
import matplotlib.pyplot as plt
ax = plt.subplot(111)
ax.set_yscale("log", nonposy='clip')
x = np.linspace(0, n_max, n_max+1)
plt.plot(x[1:],E[1:])


# In[125]:


E[1:]

