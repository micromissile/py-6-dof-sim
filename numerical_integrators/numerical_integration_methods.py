import numpy as np  

def forward_euler(f, t_s, x, h_s, vmod, amod, cmod, Auxillary_Data_Accumulated):
  """
  Performs forward Euler integration to approximate the solution of a differential equation.

  Input Args:
      f:   A function representing the right-hand side of the differential equation (dx/dt = f(t, x)).
      t_s: A vector of points in time at which numerical solutions will be approximated
      x:   The numerically approximated solution data to the DE, f
      h_s: The step size in seconds

  Returns:
      t_s: A vector of points in time at which numerical solutions was approximated
      x:   The numerically approximated solution data to the DE, f
      
  """

  for i in range(1, len(t_s)):
    dx, auxillary_data = f(t_s[i-1], x[:,i-1], vmod, amod, cmod)
    x[:,i] = x[:,i-1] + h_s * dx
    Auxillary_Data_Accumulated[:,i] = auxillary_data

  return t_s, x, Auxillary_Data_Accumulated

def AB2(f, t_s, x, h_s, vmod, amod, cmod, Auxillary_Data_Accumulated):
  """
  Performs the 2nd order Adams-Bashforth method to approximate the solution of a differential equation.

  Input Args:
      f:   A function representing the right-hand side of the differential equation (dx/dt = f(t, x)).
      t_s: A vector of points in time at which numerical solutions will be approximated
      x:   The numerically approximated solution data to the DE, f
      h_s: The step size in seconds

  Returns:
      t_s: A vector of points in time at which numerical solutions was approximated
      x:   The numerically approximated solution data to the DE, f
      
  """

  for i in range(1, len(t_s)):
    fim1, auxillary_data = f(t_s[i-1], x[:,i-1], vmod, amod, cmod)
    if i == 0:
      x[:,i] = x[:,i-1] + h_s * fim1
    else:
      fim2, _ = f(t_s[i-2], x[:,i-2], vmod, amod, cmod)
      x[:,i] = x[:,i-1] + 1.5*h_s*fim1 - 0.5*h_s*fim2
    Auxillary_Data_Accumulated[:,i] = auxillary_data

  return t_s, x, Auxillary_Data_Accumulated

def RK4(f, t_s, x, h_s, vmod, amod, cmod, Auxillary_Data_Accumulated):
  """
  Performs the 4th order Runge-Kutta method to approximate the solution of a differential equation.

  Input Args:
      f:   A function representing the right-hand side of the differential equation (dx/dt = f(t, x)).
      t_s: A vector of points in time at which numerical solutions will be approximated
      x:   The numerically approximated solution data to the DE, f
      h_s: The step size in seconds

  Returns:
      t_s: A vector of points in time at which numerical solutions was approximated
      x:   The numerically approximated solution data to the DE, f
      
  """

  for i in range(1, len(t_s)):
    fim1_k1, auxillary_data = f(t_s[i-1], x[:,i-1], vmod, amod, cmod)
    k1 = h_s*fim1_k1
    fim1_k2, _ = f(t_s[i-1] + 0.5*h_s, x[:,i-1] + 0.5*k1, vmod, amod, cmod)
    k2 = h_s*fim1_k2
    fim1_k3, _ = f(t_s[i-1] + 0.5*h_s, x[:,i-1] + 0.5*k2, vmod, amod, cmod)
    k3 = h_s*fim1_k3
    fim1_k4, _ = f(t_s[i-1] + h_s, x[:,i-1] + k3, vmod, amod, cmod)
    k4 = h_s*fim1_k4 
    x[:,i] = x[:,i-1] + 1/6*(k1 + 2.0*k2 + 2.0*k3 + k4)
    Auxillary_Data_Accumulated[:,i] = auxillary_data

  return t_s, x, Auxillary_Data_Accumulated
