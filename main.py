def vector_add(vector_a,vector_b):
  """Is the sum of the two vectors.
  
  This function takes in two vectors, representd as lists, then computes their sum and returns it as a list.

  Args:
    vector_a: An arbitrary vector of  arbitrary length. Represented by a list.
    vector_b: An arbitrary vector of the same length as vector_a. Repesented by a list.

  Returns:
    A list of the same length as vector_a which is the sum of the two input vectors.
  """

  result = [0]*len(vector_a)
  for element in range(len(vector_a)):
    result[element] = vector_a[element] + vector_b[element]
  return result
 
  
def conjugate(scalar):
  """Compute the conjugate of a complex number.

  This function takes in a complex number as its input and replaces the real part of its output with the real part of its input and replaces the imaginary part of its output with the conjuate part of the imaginary part of its input.

  Args:
    scalar: A complex number. 
  Returns:
    scalar: The complex conjugate of our input.
  """
  result = scalar.real
  result = result - scalar.imag*1j
  return result

def transpose(A):
  """Computes the Transpose of a matrix A.

  This function takes in a mxn matrix, A as its input and computes its transpose by replacing its rows to its columns and columns into rows so producing an nxm matrix, B.

  Args:
    A: a mxn matrix.
  Returns:
    B: a nxm matrix.
  """
  result = []
  for iterator in range(len(A[0]-1)):
    temp = []
    for element in range(len(A)-1):
      temp.append(A[iterator][element])
    result.append(temp)
  return result

def conjugateTranspose(A):
  """Computes the conjugate transpose of a matrix A.

  This function takes in a mxn matrix, A as its input and computes its conjugate transpose matrix as its output.

  Args:
    A: a mxn matrix.
  Returns:
    B: a nxn matrix as its conjugate transpose.
  """
  result = transpose(A)
  for iterator in range(len(A)):
   for element in  range(len(A[0]-1)):
     result[iterator][element] = conjugate(result[iterator][element])
  return result

def scalarVecMulti(scalar,vector):
  """Computes scalar vector multiplication.

  Multiplies every element of the vector by a scalar and returns the result.

  Args:
    scalar: a number..
    vector: A list of numbers respresting a vector.

  Returns:
    A list of numbers respresting a vector.
    """
  result =  [0]*len(vector)
  for iterator in range(len(result)):
    result[iterator] = scalar*vector[iterator]
  return result
 
def dot(vector01,vector02):
  """Computes the dot product.

  Takes in two vectors and computees their dot product.

  Args:
    vector01: A list of real numbers represting a vector.
    vector02: a list of real numbers of the same dimension as vector01 also respresenting a vector.

  Returns:
    the dot product of the inputs.
    """
  result = []
  for iterator in range(len(vector01)):
    result = result + conjugate(vector01[iterator]*vector02[iterator])
  return result

def two_norm(vector):
  """Calculates the 2-norm of a vector.

  Sum the squares of the elements of a given vector and returns the square root of the sum.

  Args:
    vector: A list of real numbers.
  Returns:
    A  real scalar which is the 2-norm of the given vector.
  """
  result = 0 
  for element in range(len(vector)):
    result = result + (vector[element]**2)
  result = result**(1/2)
  return result

def normalize(vector):
  """Normalize a given vector.

  Checks to see if the vector is normal or the zero vector. If not the vector is divided by its norm.

  Args:
    vector: A list of numbers respresenting a vector.

  Returns:
    A normalized vector if the input vector was not the zero vector. Prints an error otherwise.
  """

  norm = two_norm(vector)
  if (norm == 0):
    print("Invalid Input")
  elif (norm == 1):
    return vector
  else: 
    return scalarVecMulti((1/norm),vector)  


def orthognalDecomp(OrthoSet,vector):
  """Computes an orthogonal vector.

  This function computes the orthognal decomposition of vector with respect to Orthoset.

  Args:
    OrthoSet: A list of lists, where each element represents a vector and the list as       whole is orthonormal.
    vector: A vector of compatabile dimensions to the vectors in OrthoSet.

  Returns: 
    A vector which is orthogonal to the vectors in OrthoSet.
  """

  result = vector
  for iterator in range(len(OrthoSet)):
    tempScalar = dot(OrthoSet[iterator],vector)
    tempVector = scalarVecMulti(-tempScalar,OrthoSet[iterator])
    result = vectorAdd(result,tempVector)
  return result 

def van(A,x):
  """Computes the 4th degree Vandermonde matrix.

  This function computes the Vandermonde matrix by rising the vector x, to the power of the given input.

  Args:
    A: An mxn matrix.
    x: A vector.
  
  Returns:
    B: an mxn matrix.
  """
  result = []
  for exponent in range(5):
    temp = []
    for element in range(len(x)):
      temp.append(x[element]**exponent)
    result.append(temp)
  return result

def modifiedGS(A):
  """Computes the modified Gram-Schmidt matrix to get Q and R.

  This function takes in a vectors, A1,A2,..An and does the modifiedGS algotithm to compute its unitary matrix, Q and upper triangular matrix, R.

  Args:
    A: a list of n vectors.
  
  Returns:
    Q: An unitary martrix.
    R: An upper triangular matrix.
  """
  for element in range(len(A)):
    v[element] = A[element]
  for element in range(len(A)):
    r[element][element] = norm(v[element])
    q[element] = v[element]*((r[element][element])*(-1))
    for iterator in range(len(A[0]-1)):
      r[element][iterator] = dot(q[element],v[iterator])
      v[iterator] = v[iterator] + scalarVecMulti(-tempScalar,OrthoSet[iterator])
  return Q
  return R

def conjugateQ(A):
  """Computes the conjugate transpose of a matrix Q.

  This function takes in a mxn matrix, Q as its input and computes its conjugate transpose matrix as its output.

  Args:
    Q: a mxn matrix.
  Returns:
    p: a nxn matrix as its conjugate transpose.
  """
  result = transpose(Q)
  for iterator in range(len(Q)):
   for element in  range(len(Q[0]-1)):
     result[iterator][element] = conjugate(result[iterator][element])
  return result      

def sum(A,x,y):
  """Computes the summation of a series.

  This function takes in to boundary values as its inputs and a function as its inputs and computes the summation of the boundary values into the function as its output.

  Args:
    x: Lower boundary value.
    y: upper boundary value.
    A: a series.

  Returns:
  B: A series of the summation of the values.
  """
  for iterator in range(x:y)
    result[iterator] = result + result[iterator[0]-1]
    result.append(A)
  return result
  

def backSub(A,b):
  """Computes the back substitution to solve for a solution.

  This function takes in a vector,b and an upper triangular matrix, A and computes a vector,x which solves Ax=b.

  Args: 
    A: An upper triangular martrix.
    b: A vector.
  
  Returns:
    x: A vector which solves Ax=b.
  """

  result = b
  for iterator in range(len(A[0])):
    temp = (len([A[0])-1])
    element = temp-[iterator]
    result[element] = (b[element]-sum(scalarVecMulti((A[(element+1):n],result[element+1):n]))*((A[element][element])*(-1)))
  return result


def backSubR(R,p,x):
  """Computes a solution using backsub.

  This function takes in an upper triangular matrix, R and the conjugate unitary matrix,p and uses backsub to compute a solution for the system to get coefficients as solutions.

  Args:
    R: An upper triangular martrix.
    p: The of the conjugate unitary matrix.
    x: a vector

  Returns:
    H: the coefficients of the back substiution solution of the three inputs.
  """
  m = dot(p,x)
  H = backSub(R,m)

def d4Interpolent(H,x):
  """Computes the 4 degree interpolating polynomial.

  This functions takes in a variable x, and the solutions of H as the coefficients to compute the degree 4 interoplating polynomial.

  Args:
    x: The variable of the polynomial.
    H: The coefficients of the polynomial.

  Returns:
    The 4th degree interpolating polynomial.
  """
  return(H[0] + H[1]*x + H[2]*(x**2) + H[3]*(x**3) + H[4]*(x**4))