import numpy as np

def renormalize_vector(ru,rv,u,v):
  """
  Renormalize vectors to original norm

  Args:
  ru,rv : components of the transformed vector
  u,v   : components of the original vector
  dx : perturbation along x coordinate
  dy : perturbation along y coordinate

  Returns:
  ru,rv are modified in place
  """

  factor = np.sqrt ( ( u**2 + v**2 ) / ( ru**2 + rv**2 ) )
  ru[:,:] = ru[:,:] * factor
  rv[:,:] = rv[:,:] * factor

  return

def renormalize_tensor(rtxx,rtxy,rtyy,txx,txy,tyy):
  """
  Renormalize tensors to have the same eigenvalues as original tensor

  Args:
  rtxx,rtxy,rtyy : components of the transformed tensor
  txx,txy,tyy    : components of the original tensor

  Returns:
  rtxx,rtxy,rtyy are modified in place
  """

  # Apply transformation to input tensor
  #rtxx[:,:] = 
  #rtxy[:,:] = 
  #rtyy[:,:] = 
  # Still to be done

  return

