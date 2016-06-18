
################################################################################

#   Python/SymPy/GAlgebra program source code for numerical frame transformation 
# and versor decomposition into fixed-axis (Givens') rotations, together with 
# demo & test harness, for real  n-space spherical geometry  Cl(n, 0, 0) .  
# Version 1.6; date: 18/06/16; author: Fred Lunnon <Fred.Lunnon@gmail.com> 
# In command window execute: 
#   python -i /Users/fred/Desktop/frames_givens_demo/frames_givens_demo.py 

from sympy import *; 
from ga import *; 
from copy import *;  # import copy;  # fails? 
import random;  # from random import *;  # fails? 
from mv import *;  # static GA wrappers 

# Normalise versor (magnitude assumed positive) : 
#   version below fails to negotiate fl pt error? 
# def normalise (X) : return X/X.norm() # end def 
#   (symmetric) workaround below suffices for approximate versors! 
def mag2 (X) : return sympify(grade((rev(GA.mv(X))*X), 0)) # end def  # coerce! 

def normalise (X) : return X/sqrt(mag2(X)) # end def 


# Matrix transforming orthonormal frame  A  to  B  (just my little joke!), 
#   L-to-R composition, frame length = n  
def frame_transform_matrix (F, G) :  # global n; 
    return F.T*G; # end def 

# Versor transforming orthonormal frame F to G , 
#   L-to-R composition, frame length = n ; optional verbosity; 
#   optional spin continuity with previous result 
def frame_transform_versor (F, G, verb = False, Z0 = 0) : 
  # local lip, sig, sigk, disc, k, i, n, r, H, Z, R, R2, cosin, err;  
  # global GA, n, eps; 
  if verb : print; print "frame_transform_versor: F, G"; print F; print; print G; print;  # end if 
  lip = 0.7;  # angle threshold:  |lip| < 1 , here  cos pi/8 
  H = deepcopy(F); Z = GA.mv(1);  # initially  H = F , Z = 1 
  if verb : 
    print "k, sigk, cosin, R, H"; print; # end if 
  
  sig = +1;  # main loop 
  for k in range(0, n) : 
    R2 = H[k]*G[k]; cosin = sympify(grade(R2, 0));  # coerce! 
    if cosin >= lip :  # angle below threshold? 
      R = normalise(1 + R2);  # rotate 
      sigk = +1; 
    else : 
      R = normalise(G[k] - H[k]);  # reflect 
      sigk = -1; # end if 
    for i in range(k, n) :  # partially transform H 
      H[i] = sigk*grade(rev(R)*H[i]*R, 1); # end for 
    Z = Z*R; sig = sig*sigk; 
    if verb : 
      print k+1, sigk, cosin; print; print R; print; print H; print; # end if end for 
  # finally  H = (1/Z) F Z ~ +/-G , Z = prod_k R ; 
  
  # correct  sign(C) , provided  n  even! 
  if sig < 0 and n%2 == 0 : 
    Z = dual(Z); sig = +1; # end if 
  if sig < 0 : 
    print "frame_transform_versor: sign = ", sig; print; # end if 
  
  if Z0 <> 0 :  # ensure  sign(Z)  continuous 
    Z_Z = Z0*rev(Z);  # (old Z)/(new Z) ~ disc  (monomial) 
    disc = GA.mv(round(sympify(grade(Z_Z, 0))));  # coerce! 
    if disc == 0 and n%2 == 1 : 
      disc = rev(dual(GA.mv(round(sympify(grade(dual(Z_Z), 0))))));  # coerce! 
    if mag2(disc) <> 1 : 
      print "frame_transform_versor: discontinuity = ", disc; print; 
    else : Z = disc*Z; # end if end if 
  
  err = 0;  # roundoff error 
  for i in range(0, n) : 
    err = err + mag2( rev(G[i]) * (rev(Z)*F[i]*Z) - sig )/2; # end for 
  err = sqrt(err/n); 
  if err > eps : 
    print "frame_transform_versor: error = ", err; print; # end if 
  
  if verb : print "G, H, Z"; print G; print; print H; print; print Z; print; print "sig, err", sig, err; print; # end if 
  return Z; # end def 

# Convert vector to versor 
def vector_to_versor (V) : # local j, W; global n, gene; 
  W = 0; 
  for j in range(0, n) : 
    W = W + V[j]*gene[j]; # end for 
  return W; # end def 

# Convert matrix to versor: rotation angles halved; calls frame_transform_versor() ! 
#   L-to-R versors, L-to-R matrices; optional spin continuity 
def matrix_to_versor (A, Z0 = 0) : # local i; global n, gene; 
  return frame_transform_versor(gene, 
    [ vector_to_versor(A.row(i)) for i in range(0, n) ], 
    False, Z0); # end def 

# Convert versor to matrix: rotation angles doubled 
#   L-to-R versors, L-to-R matrices 
def versor_to_matrix (X) :  # local i, j; global n, gene; 
  return Matrix([ [ grade(rev(X)*gene[j]*X*gene[i], 0) 
    for j in range(0, n) ] for i in range(0, n) ]).T # end def 

# Build Givens' rotation in  ij-plane with cos(u) and sin(u) 
def givens(p, q, cu, su) :  # local i, j, R; global n; 
  R = [ [ 1-abs(sign(i-j)) 
    for j in range(0, n) ] for i in range(0, n) ];  # unit matrix 
  R[p][p] = cu; R[q][q] = cu; R[p][q] = -su; R[q][p] = su; 
  return Matrix(R); # end def 

# Factorise matrix into Givens' rotations  B = +/- R_1 ... R_m , 
#   with optional verbosity; L-to-R matrices 
def givens_factor_matrix (B, verb = False) : 
  #local A, B, C, R, i, j, k, r, s, t, err; global n, m, eps; 
  if verb : print; print "givens_factor_matrix: B"; print B; print; # end if 
  U = Matrix([ [ 1-abs(sign(i-j)) 
    for j in range(0, n) ] for i in range(0, n) ]);  # unit matrix 
  ij_k = [ [i, j] for i in range(0, n) for j in range(0, i)];  # schedule 
  
  A = B.T; C = deepcopy(U);  # initially  A = B^T , C = 1 
  R = [ U for i in range(0, m) ];  # rotation factors 
  k = -1;  # move counter 
  if verb : 
    print "k, i, j, t, s, r, R, A"; print; # end if 
  
  for j in range(0, n) : 
    for i in range(0, j) :  # main loop 
      k = k+1;  # next move 
      s = A[i, j]; t = A[i, i]; r = sqrt(s**2 + t**2); 
      R[k] = givens(i, j, t/r, s/r);  # t/r, s/r = cos(u), sin(u) 
      A = A*R[k]; C = C*R[k];  # new A_ij ~ 0 , A_jj > 0 
      if verb : 
        print k+1, i+1, j+1, t, s, r; print; print R[k]; print; print A; print; 
      # end if end for end for 
  # finally  A ~ B^T R[1] ... R[m] ~ diag(1, ..., 1, +/-1) 
  
  Z = B.T*C;  err = 0;  # roundoff error 
  for i in range(0, n) : 
    for j in range(0, n) : 
      err = err + (abs(Z[i, j]) - U[i, j])**2; # end for end for 
  err = sqrt(err/n); 
  if err > eps : print "givens_factor_matrix: error = ", err; print; # end if 
  
  if verb : print "B^T R[1] ... R[m]"; print Z; print; print "err", err; print; # end if 
  return R; # end def 

# Factorise versor into Givens' rotations: L-to-R composition; 
#   calls givens_factor_matrix();  Y = R_1 ... R_m (e_1 if odd)  
def givens_factor_versor (Y, verb = False) : 
  # local k, B, Z, M, Rmul, Rmat, sig, err; global n, m, gene, GA; 
  if verb : print; print "givens_factor_versor: Y"; print Y; print; # end if 
  if grade(Y, 0) <> 0 : M = 1;  # even grade 
  else : M = gene[0];  # odd grade 
  B = versor_to_matrix(Y*M);  
  Rmat = givens_factor_matrix(B); 
  Rmul = [ matrix_to_versor(Rmat[k]) for k in range(0, m) ] 
  
  Z = rev(Y*M);  # 1/(Y M) R_1 ... R_m  ~  1 ,  M in {1, e1} 
  for k in range(0, m) : 
    Z = Z*Rmul[k]; # end for 
  sig = grade(Z, 0); Rmul[m-1] = sig*Rmul[m-1];  # adjust spin 
  err = sqrt(mag2(Z - sig)/2); 
  if err > eps : print "givens_factor_versor: error = ", err; print; # end if 
  
  if verb : 
    print "R_1, ..., R_m"; print Rmul; print; 
    print "(1/Y) R_1 ... R_m, M"; print Z; print M; print; print "sig, err ", sig, err; print; # end if 
  return Rmul; # end def 

################################################################################

# Random real: component range (-1..+1) ; generator precision ~ 10E-16 
def rand_unit () : 
  return random.uniform(-1, +1); # end def 

# Random normalised row vector 
def rand_vector () :  # local j, ranlis, mag; global n; 
  ranlis = [ rand_unit() for j in range(0, n) ]; 
  mag = 0; 
  for j in range(0, n) : 
    mag = mag + ranlis[j]**2; # end for 
  mag = sqrt(mag); 
  return Matrix([ [ ranlis[j]/mag for j in range(0, n) ] ]); # end def 

# Random normalised versor product of  k  random vectors 
def rand_versor (k) : 
  # local R, i; global n; 
  R = 1; 
  for i in range(0, k) : 
    R = R * vector_to_versor(rand_vector()); # end for 
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

# Instantiate spherical algebra: call at top level; no default identifiers? 
gene = [];  # global generators 
n = 0; m = 0;  # vector & bivector dimension: read-only! 
def instantiate_GA (n0) :  # local n, m; global gene; 
  n = n0; m = n*(n-1)/2;  # cannot assign global  n,m ? 
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

# Verbose test suite 
def test_main () : 
  # local Amat, A, Bmat, B, Z, R, Y; global n; 
  print; print "test_main(): verbose tests"; 
  
  # Source and target frames  A,B  via random orthonormal matrices 
  Amat = rand_ortho(2*n);  # right-handed to right-handed 
  A = [ vector_to_versor(Amat.row(i)) for  i in range(0, n) ]; 
  Bmat = rand_ortho(2*n); 
  B = [ vector_to_versor(Bmat.row(i)) for  i in range(0, n) ]; 
  Z = frame_transform_versor(A, B, True); 
  Amat = rand_ortho(2*n);  # right-handed to left-handed 
  A = [ vector_to_versor(Amat.row(i)) for  i in range(0, n) ]; 
  Bmat = rand_ortho(2*n+1); 
  B = [ vector_to_versor(Bmat.row(i)) for  i in range(0, n) ]; 
  Z = frame_transform_versor(A, B, True); 
  
  # Target  B  random orthonormal matrix 
  B = rand_ortho(2*n);  # direct 
  R = givens_factor_matrix(B, True); 
  B = rand_ortho(2*n+1);  # reflective: fails for odd  n ! 
  R = givens_factor_matrix(B, True); 
  
  # Target  Y  random versor 
  Y = rand_versor(2*n);  # direct 
  R = givens_factor_versor(Y, True); 
  Y = rand_versor(2*n+1);  # reflective 
  R = givens_factor_versor(Y, True); 
  
  return "DONE"; # end def 

# Demo with  n = 2  set globally: versor spin; 
#   versor -> matrix homomorphism; spin continuity 
def spin_disc (l = 20) : 
  # local k, pifle, R, R_k, X, A, R_1, R_11, R_2, R_22; global gene; 
  (e1, e2) = gene;  # generator identifiers 
  print; print "spin_disc(): demos, l = ", l; print; 
  
  # ambiguous sign of versor, due to double cover of continuous rotation: 
  #   after  l  steps,  X_l ~ X_0  but  R_l ~ -R_0 ~ -1 ; 
  print "k, X_k, R_k"; print; 
  pifle = float(pi/l); 
  R = cos(pifle) + sin(pifle)*e1*e2;  # rotation angle 2pi/l 
  R_k = rev(R); 
  X = rand_versor(1); 
  for k in range(0, l+3) : 
    R_k = R_k*R; X = rev(R)*X*R; 
    print k; print X; print R_k; print; # end for end def 
  
  # spin continuity enforced during matrix conversion to versor; 
  #   unidirectional versor product mapping to matrix product 
  A = versor_to_matrix(R); 
  R_1 = matrix_to_versor(A, 1); # as  k = 1 
  R_11 = matrix_to_versor(A, -1); # as  k = l+1 
  R_2 = matrix_to_versor(A*A, R_1); # as  k = 2 
  R_22 = matrix_to_versor(A*A, R_11); # as  k = l+2 
  print "A, R_1, R_(l+1), R_2, R_(l+2)"; print; print A; print; print R_1; print R_11; print R_2; print R_22; print; 
  
  return "DONE"; # end def 

eps = 1.0e-12;  # rounding error bound (double-length floating-point) 
(GA, n, m) = instantiate_GA(4);  # user-selected geometry dimension 
print test_main();  # verbose tests 
(GA, n, m) = instantiate_GA(2);  # select plane geometry 
print spin_disc(); # spin demo 

################################################################################

# TODO --- 
# Wrap non-static properties: matrix.T , .row,  list.append; unit matrix, etc ?? 
# spin_disc() : demo dual continuity ?? 

