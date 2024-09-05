import numpy as np

method='batch'

def transform_tensor_2d_renorm(txx,txy,tyy,transf):
  global method

  if method == 'loop' :
    rtxx, rtxy, rtyy = transform_tensor_2d_renorm_loop(txx,txy,tyy,transf)
  elif method == 'batch' :
    rtxx, rtxy, rtyy = transform_tensor_2d_renorm_batch(txx,txy,tyy,transf)
  else:
    raise ValueError("Bad name of renormalization method")

  return rtxx, rtxy, rtyy

def transform_tensor_2d_renorm_loop(txx,txy,tyy,transf):

  # Initialize output tensors
  rtxx = np.zeros_like(txx)
  rtxy = np.zeros_like(txx)
  rtyy = np.zeros_like(txx)

  # Loop on grid points (this loop is inefficient for large problems!)
  for i in range(txx.shape[0]):
    for j in range(txx.shape[1]):
      # Extract 2D tensor
      tensor = np.array([[txx[i,j], txy[i,j]],
                        [txy[i,j], tyy[i,j]]])
      # Compute eigenvalues and eigenvectors
      eival, eivec = np.linalg.eigh(tensor)
      # Transform and renormalize first eigenvector
      ru = eivec[0,0] * transf[0,0,i,j] + eivec[1,0] * transf[0,1,i,j]
      rv = eivec[0,0] * transf[1,0,i,j] + eivec[1,0] * transf[1,1,i,j]
      norm = np.sqrt(ru**2+rv**2)
      eivec[0,0] = ru/norm
      eivec[1,0] = rv/norm
      # Transform and renormalize second eigenvector
      ru = eivec[0,1] * transf[0,0,i,j] + eivec[1,1] * transf[0,1,i,j]
      rv = eivec[0,1] * transf[1,0,i,j] + eivec[1,1] * transf[1,1,i,j]
      norm = np.sqrt(ru**2+rv**2)
      eivec[0,1] = ru/norm
      eivec[1,1] = rv/norm
      # Re-orthogonalize eigenvectors
      # dtheta = np.arctan2(u_x v_y - u_y v_x, u_x v_x + u_y v_y)

      theta1 = np.arctan2( eivec[1,0] , eivec[0,0] )
      theta2 = np.arctan2( eivec[1,1] , eivec[0,1] )
      dtheta = np.arctan2( eivec[0,0] * eivec[1,1] - eivec[1,0] * eivec[0,1] ,
                           eivec[0,0] * eivec[0,1] + eivec[1,0] * eivec[1,1] )
      ratio = eival[1] / ( eival[0] + eival[1] )
      if dtheta > 0 :
        theta1 = theta1 - ( dtheta - np.pi/2 ) * ratio
        theta2 = theta2 + ( dtheta - np.pi/2 ) * ( 1 - ratio)
      else :
        theta1 = theta1 + ( dtheta + np.pi/2 ) * ratio
        theta2 = theta2 - ( dtheta + np.pi/2 ) * ( 1 - ratio)
      delta = np.abs(theta2-theta1) - np.pi/2
      eivec[0,0] = np.cos(theta1)
      eivec[1,0] = np.sin(theta1)
      eivec[0,1] = np.cos(theta2)
      eivec[1,1] = np.sin(theta2)
      # Construct the diagonal matrix of eigenvalues
      Lambda = np.diag(eival)
      # Reconstruct the tensor with transformed eigenvectors
      tensor = eivec @ Lambda @ eivec.T
      rtxx[i,j] = tensor[0,0]
      rtxy[i,j] = tensor[0,1]
      rtyy[i,j] = tensor[1,1]

  return rtxx, rtxy, rtyy

def transform_tensor_2d_renorm_batch(txx,txy,tyy,transf):

  # Initialize output tensors
  rtxx = np.zeros_like(txx)
  rtxy = np.zeros_like(txx)
  rtyy = np.zeros_like(txx)

  # Compute eigenvalues and eigenvectors
  dtheta = np.sqrt( ( txx - tyy )**2 + 4 * txy**2 )
  eival0 = ( txx + tyy + dtheta ) / 2
  eival1 = ( txx + tyy - dtheta ) / 2
  theta0 = np.arctan2( txx - eival1, txy )
  theta1 = np.arctan2( txx - eival0, txy )
  eivec00 = np.cos(theta0)
  eivec10 = np.sin(theta0)
  eivec01 = np.cos(theta1)
  eivec11 = np.sin(theta1)

  # Transform and renormalize first eigenvector
  ru = eivec00 * transf[0,0,:,:] + eivec10 * transf[0,1,:,:]
  rv = eivec00 * transf[1,0,:,:] + eivec10 * transf[1,1,:,:]
  norm = np.sqrt(ru**2+rv**2)
  eivec00 = ru/norm
  eivec10 = rv/norm

  # Transform and renormalize second eigenvector
  ru = eivec01 * transf[0,0,:,:] + eivec11 * transf[0,1,:,:]
  rv = eivec01 * transf[1,0,:,:] + eivec11 * transf[1,1,:,:]
  norm = np.sqrt(ru**2+rv**2)
  eivec01 = ru/norm
  eivec11 = rv/norm

  # Re-orthogonalize eigenvectors
  # dtheta = np.arctan2(u_x v_y - u_y v_x, u_x v_x + u_y v_y)
  theta0 = np.arctan2( eivec10 , eivec00 )
  theta1 = np.arctan2( eivec11 , eivec01 )
  dtheta = np.arctan2( eivec00 * eivec11 - eivec10 * eivec01 ,
                       eivec00 * eivec01 + eivec10 * eivec11 )
  ratio = eival1 / ( eival0 + eival1 )
  theta0 = np.where( dtheta > 0,
                     theta0 - ( dtheta - np.pi/2 ) * ratio ,
                     theta0 + ( dtheta + np.pi/2 ) * ratio )
  theta1 = np.where( dtheta > 0,
                     theta1 + ( dtheta - np.pi/2 ) * ( 1 - ratio) ,
                     theta1 - ( dtheta + np.pi/2 ) * ( 1 - ratio) )
  delta = np.abs(theta1-theta0) - np.pi/2
  eivec00 = np.cos(theta0)
  eivec10 = np.sin(theta0)
  eivec01 = np.cos(theta1)
  eivec11 = np.sin(theta1)

  # Reconstruct the tensor with transformed eigenvectors
  rtxx = eival0 * eivec00**2 + eival1 * eivec01**2
  rtxy = eival0 * eivec00 * eivec10 + eival1 * eivec01 * eivec11
  rtyy = eival0 * eivec10**2 + eival1 * eivec11**2

  return rtxx, rtxy, rtyy

