
################################################################################

#   Python/SymPy/GAlgebra program source for numerical frame transformation 
# and versor decomposition into fixed-axis (Givens') rotations, with demo & test 
# harness, for real Clifford algebra  Cl(p,q,r) . 
# Version 2.3; date: 28/07/16; author: Fred Lunnon <Fred.Lunnon@gmail.com> 
# In command window execute: 
#   python -i /Users/fred/fred/euclidean/frames_givens_demo/frames_givens_demo.py 

from sympy import *; 
from ga import *; 
from copy import *;  # import copy;  # fails? 
import random;  # from random import *;  # fails? 
from mv import *;  # static GA wrappers 

eps = 1.0e-12;  # rounding error bound (double-length floating-point) 
FGD_version = 2.3;  # update !! 

# Normalise versor (magnitude assumed positive) : 
#   version below fails to negotiate fl pt error? 
# def normalise (X) : return X/X.norm() # end def 
#   (symmetric) workaround below suffices for approximate versors! 
def mag2 (X) : return sympify(grade((rev(GA.mv(X))*X), 0)) # end def  # coerce! 

def normalise (X) : return X/sqrt(abs(mag2(X))) # end def 


# Versor transforming orthonormal  Cl(p,q,r)  frame F to G ; 
#   L-to-R composition, frame length = n ; 
#   optional verbosity and spin continuity with previous result 
# (fails to detect mismatched isotropic  --- FIX ?? ) 
def frame_transform_versor (F, G, verb = False, Z0 = 0) : 
  # local s, t, disc, k, i, n, r, Hk, H, Z, Rk, R2, err, m1, m2, n1, n2;  
  # global GA, n, p, q, r, eps, sigs, gene; 
  if verb : 
    print; print "frame_transform_versor: F, G"; print F; print; print G; print; 
    print "k, s, ||R1||, ||R2||, R, H"; print; # end if 
  
  Z = GA.mv(1);  t = 1;  # main loop 
  for k in range(0, n) : 
    Hk = t*grade(rev(Z)*F[k]*Z, 1); 
    R1 = G[k] - Hk; m1 = mag2(R1); 
    R2 = mag2(G[k]) + Hk*G[k]; m2 = mag2(R2); 
    n1 = 0 if m1 == 0 else min(abs(m1), 1.0/abs(m1)); # end if 
    n2 = 0 if m2 == 0 else min(abs(m2), 1.0/abs(m2)); # end if 
    if max(n1, n2) < eps :  # copy (signatures clash or both isotropic) 
      Rk = GA.mv(1); s = 1; 
    elif n1 > n2 :  # reflect 
      Rk = normalise(R1); s = -sign(m1); 
    else :  # rotate 
      Rk = normalise(R2); s = +sign(m2); # end if 
    Z = Z*Rk;  t = t*s;  # accumulated transform & sign 
    if verb : 
      print k+1, s, m1, m2; print; print Rk; print; print Hk; print; # end if end for 
  # finally  H = (1/Z) F Z ~ +/-G , Z = prod_k R ; 
  
  # fix frame signs if possible ( p  even,  r  zero;  p  odd,  r  nonzero);  
  #   misses mismatched isotropic generators! 				 --- FIX ?? 
  if t < 0 : Z = Z*GA_J;  # end if 
  
  if Z0 <> 0 :  # ensure  sign(Z)  continuous 
    Z_Z = Z0*rev(Z);  # (old Z)/(new Z) ~ disc  (monomial) 
    disc = GA.mv(round(sympify(grade(Z_Z, 0))));  # coerce! 
    if disc == 0 and n%2 == 1 : 
      disc = rev(dual(GA.mv(round(sympify(grade(dual(Z_Z), 0)))))); # end if 
    if mag2(disc) <> 1 : 
      print "frame_transform_versor: discontinuity = ", disc; print; 
    else : Z = disc*Z; # end if 
  
  err = 0;  # roundoff error 
  H = [ grade(rev(Z)*F[k]*Z, 1) for k in range(0, n) ]; 
  for k in range(0, n) : 
    err = err + abs(mag2(G[k] - H[k])); # end for 
  err = sqrt(err/n).evalf(); 
  if err > eps : 
    print "frame_transform_versor: error = ", err; print; # end if 
    
  if verb : 
    print "G, H, Z"; print G; print; print H; print; print Z; print; print "t, err ", t, err; print; # end if 
  return Z; # end def 

# Convert vector to versor 
def vector_to_versor (V) :  # local j, W; global n, gene; 
  W = 0; 
  for j in range(0, n) : 
    W = W + V[j]*gene[j]; # end for 
  return W; # end def 

# Convert matrix to versor: rotation angles halved; calls frame_transform_versor() ! 
#   L-to-R versors, L-to-R matrices; optional spin continuity 
def matrix_to_versor (A, Z0 = 0) :  # local i; global n, gene; 
  return frame_transform_versor(gene, 
    [ vector_to_versor(A.row(i)) for i in range(0, n) ], 
    False, Z0); # end def 

# Convert normalised versor to matrix: rotation angles doubled; 
#   L-to-R versors, L-to-R matrices; 
#   row  i  equals vector action of versor  X  on  e_i ; 
def versor_to_matrix (X) :  # local i, j; global n, p, q, r, gene, sigs; 
  return Matrix([(rev(X)*gene[i]*X).blade_coefs(gene) for i in range(0, n)]); 
  # end def 

# Build Givens' rotation, boost, translation in  ij-plane 
def givens(i0, j0, t, s) :  # local i, j, R, r0; global n; 
  r0 = sqrt(abs(sigs[j0]*s**2 + sigs[i0]*t**2));  
  if r0 < eps+eps : r0 = 1;  # both isotropic? 
  R = [ [ 1-abs(sign(i-j)) 
    for j in range(0, n) ] for i in range(0, n) ];  # unit matrix 
  R[i0][i0] = t/r0; R[j0][j0] = t/r0;  # cos(angle)  etc. 
  R[i0][j0] = -s/r0; R[j0][i0] = sigs[i0]*sigs[j0]*s/r0;  # sin(angle)  etc. 
  return Matrix(R); # end def 

# Factorise O(p,q,r) matrix into Givens' rotations  B = +/- R_1 ... R_m , 
#   with optional verbosity; L-to-R matrices 
def givens_factor_matrix (B, verb = False) : 
  #local A, C, R, i, j, k, r, s, t, err; global n, m, eps; 
  if verb : print; print "givens_factor_matrix: B"; print B; print; # end if 
  
  C = Matrix([ [ 1-abs(sign(i-j)) 
    for j in range(0, n) ] for i in range(0, n) ]);  # unit matrix 
  A = deepcopy(B); # initially  A = B , C = 1 
  R = [ C for i in range(0, m) ];  # rotation factors 
  k = m;  # move counter 
  if verb : 
    print "k, i, j, R, A"; print; # end if 
  
  for j in range(0, n) : 
    for i in range(0, j) : 
      k = k-1;  # main loop: next move 
      R[k] = givens(i, j, A[i, i], A[i, j]); 
      A = A*R[k];  # new A_ij ~ 0 , A_jj > 0 
      R[k][i, j] = -R[k][i, j];  R[k][j, i] = -R[k][j, i]; 
      C = R[k]*C;  # update product with inverse G-rot 
      if verb : 
        print k+1, i+1, j+1; print; print R[k]; print; print A; print; 
      # end if end for end for 
  # finally  A ~ diag(1, ..., 1, +/-1) , C ~ B (except row  n ?) 
  
  err = 0;  # roundoff error 
  for i in range(0, n) : 
    for j in range(0, n) : 
      err = err + (B[i, j] - C[i, j])**2; # end for end for 
  err = sqrt(err/n).evalf();  # unevaluated? 
  if err > eps : print "givens_factor_matrix: error = ", err; print; # end if 
  
  if verb : print "B, R[1] ... R[m]"; print B; print; print C; print; print "err", err; print; # end if 
  return R; # end def 

# Factorise  Cl(p,q,r)  versor into Givens' rotations: L-to-R composition; 
#   calls givens_factor_matrix();  Y = R_1 ... R_m (e_1 if odd) 
def givens_factor_versor (Y, verb = False) : 
  # local k, B, Z, M, Rmul, Rmat, sig, err; global n, m, gene, GA; 
  if verb : print; print "givens_factor_versor: Y"; print Y; print; # end if 
  if grade(Y, 0) <> 0 : M = 1;  # even grade 
  else : M = gene[0];  # odd grade:  e_1  non-isotropic? 
  B = versor_to_matrix(Y*M); 
  Rmat = givens_factor_matrix(B); 
  Rmul = [ matrix_to_versor(Rmat[k]) for k in range(0, m) ] 
  
  Z = rev(Y*M);  # 1/(Y M) R_1 ... R_m  ~  1 ,  M in {1, e_1} 
  for k in range(0, m) : 
    Z = Z*Rmul[k]; # end for 
  sig = grade(Z, 0); Rmul[m-1] = sig*Rmul[m-1];  # adjust spin 
  err = sqrt(abs(mag2(Z - sig))/2).evalf(); 
  if err > eps : print "givens_factor_versor: error = ", err; print; # end if 
  
  if verb : 
    print "R_1, ..., R_m"; print Rmul; print; 
    print "(1/Y) R_1 ... R_m, M"; print Z; print M; print; print "sig, err ", sig, err; print; # end if 
  return Rmul; # end def 

################################################################################

# Random real: component range (-1..+1) ; generator precision ~ 10E-16 
def rand_unit () : 
  return random.uniform(-1, +1); # end def 

# Random real: component range (-1..+1) ; generator precision ~ 10E-16 
def rand_unit () : 
  return random.uniform(-1, +1); # end def 

# Random row vector (non-isotropic) 
def rand_vector () :  # local j, ranlis, mag; global n, sigs; 
  ranlis = [ rand_unit() for j in range(0, n) ]; 
  return Matrix([ [ ranlis[j] for j in range(0, n) ] ]); # end def 

# Random orthonormal  Cl(p,q,r)  versor product of  l  vectors, avoiding near-isotropic 
def rand_versor (l) : 
  # local X, R, i; global n; 
  if l%2 == 0 : X = 1 
  else : X = vector_to_versor(rand_vector()); # end if 
  for i in range(0, floor(l/2)) : 
    X = X * (1 + vector_to_versor(rand_vector()) * 
      vector_to_versor(rand_vector())); # end for 
  return normalise(X); # end def 

# Random orthonormal  O(p,q,r)  matrix product of  l  reflections 
#   Q & D : for  n  odd  &  l  odd, version below yields  det = +1 ? 
def rand_ortho (l) : 
  return versor_to_matrix(rand_versor(l)); # end def 

# Instantiate Clifford algebra, given signature list comprising  +1,-1,0's 
def instantiate_GA (sigs0) :  # local j; 
  global GA, GA_J, gene, sigs, eps, n, m; 
  global e1, e2, e3, e4, e5, e6;  # export default identifiers 
  sigs = sigs0; 
  n = len(sigs0); m = n*(n-1)/2; 
  
  if n == 2 : 
    GA = Ga('e1 e2', g = sigs); 
    (e1, e2) = GA.mv(); gene = [e1, e2]; 
  elif n == 3 : 
    GA = Ga('e1 e2 e3', g = sigs); 
    (e1, e2, e3) = GA.mv(); gene = [e1, e2, e3]; 
  elif n == 4 : 
    GA = Ga('e1 e2 e3 e4', g = sigs); 
    (e1, e2, e3, e4) = GA.mv(); gene = [e1, e2, e3, e4]; 
  elif n == 5 : 
    GA = Ga('e1 e2 e3 e4 e5', g = sigs); 
    (e1, e2, e3, e4, e5) = GA.mv(); gene = [e1, e2, e3, e4, e5]; 
  elif n == 6 : 
    GA = Ga('e1 e2 e3 e4 e5 e6', g = sigs); 
    (e1, e2, e3, e4, e5, e6) = GA.mv(); gene = [e1, e2, e3, e4, e5, e6]; 
  else : 
    print; print "You're on your own, sunshine!  n = ", n; print; # end if 
  
  GA_J = 1;  # quasi-pseudar 
  for j in range(0, n) : 
    if sigs[j] <> 0 : GA_J = GA_J*gene[j]; # end if end for 
  return None; # end def 

# Verbose test suite 
def test_main () : 
  # local Amat, A, Bmat, B, Z, R, Y; global n, m, p, q, r; 
  print; print "test_main(): verbose tests"; 
  
  # Source and target frames  A,B  via random orthonormal matrices 
  Amat = rand_ortho(floor(n/2)*2);  # right-handed to right-handed 
  A = [ vector_to_versor(Amat.row(i)) for i in range(0, n) ]; 
  Bmat = rand_ortho(floor(n/2)*2); 
  B = [ vector_to_versor(Bmat.row(i)) for i in range(0, n) ]; 
  Z = frame_transform_versor(A, B, True); 
  Amat = rand_ortho(floor(n/2)*2);  # right-handed to left-handed 
  A = [ vector_to_versor(Amat.row(i)) for i in range(0, n) ]; 
  Bmat = rand_ortho(floor(n/2)*2+1); 
  B = [ vector_to_versor(Bmat.row(i)) for i in range(0, n) ]; 
  Z = frame_transform_versor(A, B, True); 
  
  # Target  B  random orthonormal matrix 
  B = rand_ortho(floor(n/2)*2);  # direct 
  R = givens_factor_matrix(B, True); 
  B = rand_ortho(floor(n/2)*2+1);  # reflective: fails for odd  n ! 
  R = givens_factor_matrix(B, True); 
  
  # Target  Y  random (normalised?) versor 
  Y = rand_versor(floor(n/2)*2);  # direct 
  R = givens_factor_versor(Y, True); 
  Y = rand_versor(floor(n/2)*2+1);  # reflective 
  R = givens_factor_versor(Y, True); 
  
  return "DONE"; # end def 

# Demo with  n = 2  set globally: versor spin; 
#   versor -> matrix homomorphism; spin continuity 
def spin_disc (l = 20) : 
  # local k, pile, R, R_k, X, A, R_1, R_11, R_2, R_22; global gene; 
  (e1, e2) = gene;  # generator identifiers 
  print; print "spin_disc(): demos, l = ", l; print; 
  
  # ambiguous sign of versor, due to double cover of continuous rotation: 
  #   after  l  steps,  X_l ~ X_0  but  R_l ~ -R_0 ~ -1 ; 
  print "k, X_k, R_k"; print; 
  pile = float(pi/l); 
  R = cos(pile) + sin(pile)*e1*e2;  # rotation angle 2pi/l 
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

# Time 02:21 
print; print "Python/SymPy/GAlgebra frames_givens_demo version", FGD_version; 
instantiate_GA([1,1,1,1]);  # 3-sphere geometry  Cl(4) 
print test_main();  # verbose tests 
instantiate_GA([1,1]);  # circle geometry  Cl(2) 
print spin_disc(); # spin demo 
instantiate_GA([1,1,-1,-1,0,0]);  # mixed degenerate  Cl(2,2,2) 
print test_main();  # verbose tests 

################################################################################

# TODO --- 
# Wrap non-static properties: matrix.T , .row, list.append; unit matrix, etc ?? 
# spin_disc() : demo dual continuity ?? 
# test_main () : tidy identifiers ?? 
# rand_versor() : fix sign quibble? 
# frame_transform_versor() : 
#   misses mismatched isotropic generators! 
# givens_factor_matrix() : 
#   query: test translation in 1-space CGA = GA(2, 1) ? 
#   Cl(2,2,2) second call --- error =  0.816496580927725 ?? 
# Version 3 : arithmetic wrappers throughout; latter GA_multor, GAlgebra, Clifford ... 
