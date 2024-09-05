"""
Tools to transform vectors and tensors following the perturbation of the coordinates

Available functions:
 - transform_vector: transform vector to follow the perturbation of the coordinates
 - transform_tensor: transform symmetric tensor to follow the perturbation of the coordinates

Module parameters:

"""
import numpy as np
import util_grid

# Method used for transformation:
# jacobian : apply full jacobian transformation
# renormalization : apply full jacobian transformation + renormalization
# rotation_only : apply only rotation component of jacobian
method = 'renormalization'

def transform_vector(u,v,dx,dy,grid_type='T'):
  """
  Transform vectors to follow the perturbation of the coordinates

  Args:
  u,v : components of the vector
  dx : perturbation along x coordinate
  dy : perturbation along y coordinate

  Returns:
  u,v are modified in place
  """

  # Dimension of input variable
  ndim = u.ndim

  if grid_type == 'UV' :
    # Calculate transformation at U and V points
    transf_u = compute_transformation(dx,dy,grid_type='U')
    transf_v = compute_transformation(dx,dy,grid_type='V')
  else:
    # Calculate transformation at required grid location
    transf_u = compute_transformation(dx,dy,grid_type)
    transf_v = compute_transformation(dx,dy,grid_type)

  # Transform vector for every 2d slice
  if ndim == 2 :
    u2d = u[:,:]
    v2d = v[:,:]
    ru2d, rv2d = transform_vector_2d(u2d,v2d,transf_u,transf_v,grid_type)
    u[:,:] = ru2d
    v[:,:] = rv2d
  elif ndim == 3 :
    for k in range(u.shape[0]):
      u2d = u[k,:,:]
      v2d = v[k,:,:]
      ru2d, rv2d = transform_vector_2d(u2d,v2d,transf_u,transf_v,grid_type)
      u[k,:,:] = ru2d
      v[k,:,:] = rv2d
  elif ndim == 4 :
    for l in range(u.shape[0]):
      for k in range(u.shape[1]):
        u2d = u[l,k,:,:]
        v2d = v[l,k,:,:]
        ru2d, rv2d = transform_vector_2d(u2d,v2d,transf_u,transf_v,grid_type)
        u[l,k,:,:] = ru2d
        v[l,k,:,:] = rv2d
  else:
    raise ValueError("Bad variable dimension")

def transform_vector_2d(u,v,transf_u,transf_v,grid_type='T'):
  """
  Transform vectors to follow the perturbation of the coordinates

  Args:
  u,v : components of the vector
  transf_u, transf_v : transformation at U and V points

  Returns:
  ru,rv : transformed components of vector
  """

  # Define working arrays
  u_v = np.zeros_like(u)
  v_u = np.zeros_like(v)
  ru  = np.zeros_like(u)
  rv  = np.zeros_like(v)

  if grid_type == 'UV' :
    # Compute u at V points
    u_v[1:,1:] = ( u[0:-1, 0:-1] + u[1:, 0:-1] + u[0:-1, 1:] + u[1:, 1:] ) / 4
    u_v[0,:] = u[0,:] ; u_v[:,0] = u[:,0]
    # Compute v at U points
    v_u[1:,1:] = ( v[0:-1, 0:-1] + v[1:, 0:-1] + v[0:-1, 1:] + v[1:, 1:] ) / 4
    v_u[0,:] = v[0,:] ; v_u[:,0] = v[:,0]
    # Apply transformation to input vector
    ru = u   * transf_u[0,0,:,:] + v_u * transf_u[0,1,:,:]
    rv = u_v * transf_v[1,0,:,:] + v   * transf_v[1,1,:,:]
    # Renormalize vector if required by user
    if method == 'renormalization' :
      rv_u = u   * transf_v[1,0,:,:] + v_u * transf_v[1,1,:,:]
      ru_v = u_v * transf_u[0,0,:,:] + v   * transf_u[0,1,:,:]
      renormalize_vector(ru  ,rv_u,u  ,v_u)
      renormalize_vector(ru_v,rv  ,u_v,v  )
  else:
    # Apply transformation to input vector
    ru = u * transf_u[0,0,:,:] + v * transf_u[0,1,:,:]
    rv = u * transf_v[1,0,:,:] + v * transf_v[1,1,:,:]
    # Renormalize vector if required by user
    if method == 'renormalization' :
      renormalize_vector(ru,rv,u,v)

  return ru, rv

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

def transform_tensor(txx,txy,tyy,dx,dy,grid_type='T'):
  """
  Transform symmetric tensor to follow the perturbation of the coordinates

  Args:
  txx,txy,tyy : components of the tensor
  dx : perturbation along x coordinate
  dy : perturbation along y coordinate
  grid_type : type of staggered grid (T, U, V or F)

  Returns:
  rtxx,rtxy,rtyy : transformed components of the tensor
  """

  # Dimension of input variable
  ndim = txx.ndim

  # Calculate the transformation at appropriate grid points
  transf = compute_transformation(dx,dy,grid_type)

  # Transform tensor for every 2d slice
  if ndim == 2 :
    txx2d = txx[:,:]
    txy2d = txy[:,:]
    tyy2d = tyy[:,:]
    rtxx2d, rtxy2d, rtyy2d = transform_tensor_2d(txx2d,txy2d,tyy2d,transf)
    txx[:,:] = rtxx2d
    txy[:,:] = rtxy2d
    tyy[:,:] = rtyy2d
  elif ndim == 3 :
    for k in range(txx.shape[0]):
      txx2d = txx[k,:,:]
      txy2d = txy[k,:,:]
      tyy2d = tyy[k,:,:]
      rtxx2d, rtxy2d, rtyy2d = transform_tensor_2d(txx2d,txy2d,tyy2d,transf)
      txx[k,:,:] = rtxx2d
      txy[k,:,:] = rtxy2d
      tyy[k,:,:] = rtyy2d
  elif ndim == 4 :
    for l in range(txx.shape[0]):
      for k in range(txx.shape[1]):
        txx2d = txx[l,k,:,:]
        txy2d = txy[l,k,:,:]
        tyy2d = tyy[l,k,:,:]
        rtxx2d, rtxy2d, rtyy2d = transform_tensor_2d(txx2d,txy2d,tyy2d,transf)
        txx[l,k,:,:] = rtxx2d
        txy[l,k,:,:] = rtxy2d
        tyy[l,k,:,:] = rtyy2d
  else:
    raise ValueError("Bad variable dimension")

  return

def transform_tensor_2d(txx,txy,tyy,transf):
  """
  Transform symmetric tensor to follow the perturbation of the coordinates

  Args:
  txx,txy,tyy : components of the tensor
  transf : transformation to apply

  Returns:
  rtxx,rtxy,rtyy : transformed components of the tensor
  """

  if method == 'renormalization' :
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
  else:
    # Directly apply transformation to input tensor
    rtxx = txx * transf[0,0,:,:]**2 + 2 * txy * transf[0,0,:,:] * transf[0,1,:,:] \
           + tyy * transf[0,1,:,:]**2
    rtxy = txx * transf[1,0,:,:] * transf[0,0,:,:] + tyy * transf[0,1,:,:] * transf[1,1,:,:] \
           + txy * ( transf[0,1,:,:] * transf[1,0,:,:] + transf[0,0,:,:] * transf[1,1,:,:] )
    rtyy = txx * transf[1,0,:,:]**2 + 2 * txy * transf[1,0,:,:] * transf[1,1,:,:] \
           + tyy * transf[1,1,:,:]**2

  return rtxx, rtxy, rtyy

def compute_transformation(dx,dy,grid_type='T'):
  """
  Compute transformation associated to perturbation

  Args:
  dx : perturbation along x coordinate
  dy : perturbation along y coordinate
  grid_type : type of staggered grid (T, U, V or F)

  Returns:
  transf : transformation
  """

  # Initialize transformation
  transf = np.zeros( (2, 2) + dx.shape )

  # Set shifts to use for each type of grid
  if grid_type == 'T':
    ish = 2 ; jsh = 2
  elif grid_type == 'U':
    ish = 1 ; jsh = 2
  elif grid_type == 'V':
    ish = 2 ; jsh = 1
  elif grid_type == 'F':
    ish = 1 ; jsh = 1
  else:
    raise ValueError("Bad type of staggered grid (T, U , V or F)")

  # We remove translation by computing difference of dx and dy
  # (dx and dy are always assumed to correspond to grid T)
  ilo = ish - 1 ; jlo = jsh - 1
  dx_x = ( dx[jlo:-1,ish:] - dx[jlo:-1,:-ish] ) / ish
  dx_y = ( dx[jsh:,ilo:-1] - dx[:-jsh,ilo:-1] ) / jsh
  dy_x = ( dy[jlo:-1,ish:] - dy[jlo:-1,:-ish] ) / ish
  dy_y = ( dy[jsh:,ilo:-1] - dy[:-jsh,ilo:-1] ) / jsh

  # The Jacobian of the transformation of the coordinates is then J = I + dx'/dx
  if method == 'rotation_only' :
    # J is the combination of a rotation and deformation: J = R D.
    # The rotation angle is computed to obtain a symmetric matrix D from T:
    # tan(transf) = [ (Tyx - Txy)/2 ] / [ 1 + (Txx + Tyy)/2 ]
    trace = ( dx_x + dy_y ) / 2
    asymm = ( dy_x - dx_y ) / 2
    theta = np.arctan2( asymm, 1 + trace )
    # The transformation is then a rotation by theta
    transf[0,0,jlo:-1,ilo:-1] = np.cos(theta)
    transf[0,1,jlo:-1,ilo:-1] = np.sin(theta)
    transf[1,0,jlo:-1,ilo:-1] = - np.sin(theta)
    transf[1,1,jlo:-1,ilo:-1] = np.cos(theta)
  else:
    # We apply the inverse Jacobian transformation
    det = ( 1 + dx_x ) * ( 1 + dy_y ) - dx_y * dy_x
    transf[0,0,jlo:-1,ilo:-1] = ( 1 + dy_y ) / det
    transf[0,1,jlo:-1,ilo:-1] = - dx_y / det
    transf[1,0,jlo:-1,ilo:-1] = - dy_x / det
    transf[1,1,jlo:-1,ilo:-1] = ( 1 + dx_x ) / det
    # Check definite positiveness of the Jacobian
    count = np.count_nonzero( (1+dx_x > 0) & (det > 0) )
    if count < det.size :
      formatted_count = f"{100*(1-count/det.size):.2f}"
      print("E R R O R :")
      print(f"Transformation is singular on {formatted_count}% of the domain")
      print("Consider reducing the ratio between the standard deviation and the correlation length scale")
      print("or using the rotation_only option for vectors and tensors.")
      raise ValueError("Singular transformation") 

  # Extrapolate to boundaries
  for i in range(1):
    for j in range(1):
      transf[i,j,:ilo,:] = transf[i,j,ilo,:][np.newaxis, :]
      transf[i,j,-1,:]   = transf[i,j,-2,:]
      transf[i,j,:,:jlo] = transf[i,j,:,jlo][:, np.newaxis]
      transf[i,j,:,-1]   = transf[i,j,:,-2]

  return transf

