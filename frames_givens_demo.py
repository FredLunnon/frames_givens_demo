
################################################################################

# Frame conversion and factoring into Givens' rotations for versors; 
#   left-to-right composition. 
# Version 1.2; date: 12/06/16; author: Fred Lunnon <Fred.Lunnon@gmail.com> 
# In command window execute: 
#   python -i /Users/fred/fred/euclidean/frames_givens_demo.py 

from sympy import *; 
from ga import *; 
from copy import *;  # import copy;  # fails? 
import random;  # from random import *;  # fails? 

# Normalise versor (magnitude assumed positive) : 
#   version below fails to negotiate fl pt error? 
# def normalise (X) : return X/X.norm() # end def 
#   (symmetric) workaround below suffices for approximate versors! 
def mag2 (X) : return sympify((X.rev()*X).grade(0)) # end def  # coerce! 

def normalise (X) : return X/sqrt(mag2(X)) # end def 

# Versor transforming orthonormal frame A to B , 
#   with optional verbosity; L-to-R composition, frame length = n 
def frame_transform (A, B, verb = False) : 
  # local lip, sig, sigk, k, i, n, r, C, Z, R, R2, cosin; # global GA, n, eps; 
  if verb : print; print "frame_transform: A, B"; print A; print; print B; print;  # end if 
  lip = 0.7;  # angle borderline: |lip| < 1 , here  cos pi/8 
  C = deepcopy(A); Z = GA.mv(1);  # initially  C = A , Z = 1 
  if verb : 
    print "k, sigk, cosin, R, C"; print; # end if 
  
  sig = +1;  # main loop  
  for k in range(0, n) : 
    R2 = C[k]*B[k]; cosin = sympify(R2.grade(0));  # coerce! 
    if cosin > lip :  # inspect angle 
      R = normalise(1 + R2);  # rotate 
      sigk = +1; 
    else : 
      R = normalise(B[k] - C[k]);  # reflect 
      sigk = -1;  # end if 
    for i in range(k, n) :  # partially transform C 
      C[i] = sigk*(R.rev()*C[i]*R).grade(1); # end for 
    Z = Z*R; sig = sig*sigk; 
    if verb : 
      print k+1, sigk, cosin; print; print R; print; print C; print; # end if end for 
  
  # finally  C = (1/Z) A Z = +/-B , Z = prod_k R ; fix sign ( n  even) 
  if sig < 0 and n%2 == 0 : 
    Z = Z.dual(); sig = +1; # end if 
  err = 0;  # inspect roundoff error 
  for i in range(0, n) : 
    err = err + mag2( B[i].rev() * (Z.rev()*A[i]*Z) - sig )/2;  # end for  
  err = sqrt(err/n); 
  if err > eps : 
    print "frame_transform: error = ", err; print; 
  elif sig < 0 :  
    print "frame_transform: sign = ", sig; print; # end if 
  
  if verb : print "B, C, Z"; print B; print; print C; print; print Z; print; print "sig, err", sig, err; print; # end if 
  return Z; # end def 

# Convert vector to versor 
def vector_to_versor (V) : # local j, W; global n, gene; 
  W = 0; 
  for j in range(0, n) : 
    W = W + V[j]*gene[j];  # end for 
  return W; # end def 

# Convert matrix to versor: rotation angles halved; calls frame_transform() ! 
#   L-to-R versors, L-to-R matrices 
def matrix_to_versor (A) : # local i; global n, gene; 
  return frame_transform(gene, 
    [ vector_to_versor(A.row(i)) for i in range(0, n) ]); # end def 

# Convert versor to matrix: rotation angles doubled 
#   L-to-R versors, L-to-R matrices 
def versor_to_matrix (X) :  # local i, j; global n, gene; 
  return Matrix([ [ (X.rev()*gene[j]*X*gene[i]).grade(0) 
    for j in range(0, n) ] for i in range(0, n) ]).T # end def 

# Build Givens' rotation in  ij-plane with cos(u) and sin(u) 
def given(p, q, cu, su) :  # local i, j, R; global n; 
  R = [ [ 1-abs(sign(i-j)) 
    for j in range(0, n) ] for i in range(0, n) ];  # unit matrix 
  R[p][p] = cu; R[q][q] = cu; R[p][q] = -su; R[q][p] = su; 
  return Matrix(R); # end def 

# Factorise matrix into Givens' rotations  B = R_1 ... R_m , 
#   with optional verbosity; L-to-R matrices 
def givens_factor_matrix (B, verb = False) : 
  #local A, B, C, R, i, j, k, r, s, t, err; global n, m, eps; 
  if verb : print; print "givens_factor_matrix: B"; print B; print;  # end if 
  U = Matrix([ [ 1-abs(sign(i-j)) 
    for j in range(0, n) ] for i in range(0, n) ]);  # unit matrix 
  ij_k = [ [i, j] for i in range(0, n) for j in range(0, i)]; # schedule 
  
  A = B.T; C = deepcopy(U);  # initially  A = B^T , C = 1 
  R = [ U for i in range(0, m) ];  # rotation factors 
  if verb : 
    print "k, i, j, t, s, r, R, A"; print;  # end if 
  
  for k in range(0, m) :  # main loop 
    i = ij_k[k][0]; j = ij_k[k][1];  # next move 
    # Givens' rotation clearing matrix coefficient  [i, j]  from  A : 
    #   t/r, s/r = cos u, sin u ; 
    s = A[j, i]; t = A[j, j]; r = sqrt(s**2 + t**2); 
    R[k] = given(j, i, t/r, s/r); # t,s = cos,sin(angle) 
    A = A*R[k]; C = C*R[k]; # new A_ij ~ 0 
    if verb : 
      print k+1, i+1, j+1, t, s, r; print; print R[k]; print; print A; print; # end if end for 
  Z = B.T*C  # finally  A ~ diag(1, ..., 1, +/-1) ,  C ~ B  (if special) 
  
  # B^T R[1] ... R[m] ~ I, diag(1, ..., 1, -1)  for special, reflecting resp. 
  err = 0;  
  for i in range(0, n) : 
    for j in range(0, n) : 
      err = err + (abs(Z[i, j]) - U[i, j])**2;  # end for end for 
  err = sqrt(err/n); 
  if err > eps : print "givens_factor: error = ", err; print;  # end if 

  if verb : print "B^T R[1] ... R[m]"; print Z; print; print "err", err; print;  # end if 
  return R; # end def 

# Factorise versor into Givens' rotations: L-to-R composition; 
#   calls givens_factor_matrix();  Y = R_1 ... R_m M 
def givens_factor_versor (Y, verb = False) :  
  # local k, B, Z, M, Rmul, Rmat; global n, m, gene, GA; 
  if verb : print; print "givens_factor_versor: Y"; print Y; print;  # end if 
  B = versor_to_matrix(Y); 
  Rmat = givens_factor_matrix(B); 
  Rmul = [ matrix_to_versor(Rmat[k]) for k in range(0, m) ]  
  
  Z = Y.rev();  # (1/Y) R_1 ... R_m  ~  +/-M , monomial 
  for k in range(0, m) : 
    Z = Z*Rmul[k];  # end for  
  if Z.grade(0) <> 0 : M = GA.mv(1).rev();  # expected monomial: grade even 
  elif n%2 <> 0 : M = GA.mv(1).dual().rev();  # grade odd, n odd 
  else : M = gene[n-1].dual().rev();  # end if  # grade odd, n even 
  sig = (M*Z).grade(0); # correct for spin continuity 
  Rmul[m-1] = sig*Rmul[m-1];  
  err = sqrt(mag2(sig*Z - M.rev())/2); 
  if err > eps : print "givens_factor: error = ", err; print;  # end if 
  if M <> 1 : print "givens_factor: monomial = ", M; print;  # end if 
  
  if verb : 
    print "R_1, ..., R_m"; print Rmul; print; 
    print "(1/Y) R_1 ... R_m, M"; print Z; print M; print; print "sig, err ", sig, err; print;  # end if 
  return Rmul; # end def 

################################################################################

# Random real, length-n  vector, grade-k  versor, orthonormal matrix; 
#   component range (-1..+1) ; generator reciprocal precision ~ 10^16 
def rand_unit () : 
  return random.uniform(-1, +1); # end def 

# random normalised vector qua row matrix 
def rand_vector () :  # local j, ranlis, mag; global n; 
  ranlis = [ rand_unit() for j in range(0, n) ]; 
  mag = 0; 
  for j in range(0, n) : 
    mag = mag + ranlis[j]**2;  # end for 
  mag = sqrt(mag); 
  return Matrix([ [ ranlis[j]/mag for j in range(0, n) ] ]); # end def 

# random normalised versor product of  k  random vectors 
def rand_versor (k) : 
  # local R, i; global n; 
  R = 1; 
  for i in range(0, k) : 
    R = R * vector_to_versor(rand_vector());  # end for 
  return R; # end def 

# Random orthonormal matrix product of  l  reflections 
#   For  n  odd  & product odd, version below yields  det = 1 ? 
# def rand_ortho (l) : 
#   return versor_to_matrix(normalise(rand_versor(l))); # end def 
#   purely matrix-based alternative below!  
def rand_ortho (l) : 
  # local k, R, F, A, U; global n; 
  U = Matrix([ [ 1-abs(sign(i-j)) 
    for j in range(0, n) ] for i in range(0, n) ]);  # unit matrix 
  A = U; 
  for k in range(0, l) : 
    F = rand_vector(); # row 
    R = U - 2*F.T*F;  # R = I - 2 F'*F 
    A = R*A;  # end for 
  return A; # end def 

# Instantiate spherical algebra --- no default identifiers? 
gene = [];  # global generators 
n = 0; m = 0; # vector & bivector dimension: read-only! 
def instantiate_GA (n0) : # local n, m; global gene; 
  n = n0; m = n*(n-1)/2; # cannot assign global  n,m ? 
  while len(gene) > 0 : gene.pop();  # delete generators 
  
  if n == 2 : 
    GA = Ga('e1 e2', g = [1 for i in range(0, n)]); 
    (e1, e2) = GA.mv(); 
    gene.append(e1); gene.append(e2);  # etc. below ?? 
  elif n == 3 : 
    GA = Ga('e1 e2 e3', g = [1 for i in range(0, n)]); 
    (e1, e2, e3) = GA.mv(); 
    gene.append(e1); gene.append(e2); gene.append(e3); 
  elif n == 4 : 
    GA = Ga('e1 e2 e3 e4', g = [1 for i in range(0, n)]); 
    (e1, e2, e3, e4) = GA.mv(); 
    gene.append(e1); gene.append(e2); gene.append(e3); gene.append(e4); 
  elif n == 5 : 
    GA = Ga('e1 e2 e3 e4 e5', g = [1 for i in range(0, n)]); 
    (e1, e2, e3, e4, e5) = GA.mv(); 
    gene.append(e1); gene.append(e2); gene.append(e3); gene.append(e4); gene.append(e5); 
  elif n == 6 : 
    GA = Ga('e1 e2 e3 e4 e5 e6', g = [1 for i in range(0, n)]); 
    (e1, e2, e3, e4, e5, e6) = GA.mv(); 
    gene.append(e1); gene.append(e2); gene.append(e3); gene.append(e4); gene.append(e5); gene.append(e6); 
  else : 
    print; print "You're on your own, sunshine!  n = ", n; print; # end if 
  return (GA, n, m); # end def 

# Test suite 
def test_main () : 
  
  # Source and target frames  A,B  via random orthonormal matrices 
  Amat = rand_ortho(2*n);  # right-handed to right-handed 
  A = [ vector_to_versor(Amat.row(i)) for  i in range(0, n) ]; 
  Bmat = rand_ortho(2*n); 
  B = [ vector_to_versor(Bmat.row(i)) for  i in range(0, n) ]; 
  Z = frame_transform(A, B, True); 
  Amat = rand_ortho(2*n);  # right-handed to left-handed 
  A = [ vector_to_versor(Amat.row(i)) for  i in range(0, n) ]; 
  Bmat = rand_ortho(2*n+1); 
  B = [ vector_to_versor(Bmat.row(i)) for  i in range(0, n) ]; 
  Z = frame_transform(A, B, True);  # fails: error = sqrt2 
  
  # Target  B  random orthonormal matrix 
  B = rand_ortho(2*n);  # direct 
  R = givens_factor_matrix(B, True); 
  B = rand_ortho(2*n+1);  # reflective 
  R = givens_factor_matrix(B, True); 
  
  # Target  Y  random versor 
  Y = rand_versor(2*n); # direct 
  R = givens_factor_versor(Y, True); 
  Y = rand_versor(2*n+1);  # reflective 
  R = givens_factor_versor(Y, True); 
  
  return "DONE"; # end def 

eps = 1.0e-12; # acceptable error using double-length floating-point 
(GA, n, m) = instantiate_GA(4);  # user select geometry 
print(test_main());  # run demos 

################################################################################

# TODO --- 
# Detailed specs !! 
# frame_transform () : ensure continuity of  sign(Z)  ?? 
# test_main () : demo versor -> matrix homomorphism; demo spin ?? 
# Wrap non-static properties: matrix.T , multor.rev() , unit matrix etc ?? 
# Tidy identifiers: matrix  Xm  etc, versor  Xv , frame  Xs , etc? 	?? 
#   Or  matrix & vector  A, B, C, ... ; versor & vector  ..., X, Y, Z . 

