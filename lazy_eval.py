
################################################################################

# Lazy evaluation test, using  frames_givens_demo.py 

if true : 
  (GA, n, p, q, r, m) = instantiate_GA(4, 0, 0);  # 3-sphere geometry 
  B = rand_ortho(floor(n/2)*2);  
  C = rand_ortho(floor(n/2)*2); 
  
  err = 0; 
  for i in range(0, n) : 
    for j in range(0, n) : 
      err = err + (B[i, j] - C[i, j])**2; 
  err; # unevaluated? 
  err.evalf(); 
  
  err = 0; 
  for i in range(0, n) : 
    for j in range(0, n) : 
      err = err + B[i, j]**2; 
  err; # evaluated! 
  
  B = Matrix([ [ sympify(B[i, j]) 
    for j in range(0, n) ] for i in range(0, n) ]); 
  C = Matrix([ [ sympify(C[i, j]) 
    for j in range(0, n) ] for i in range(0, n) ]); 
  
  err = 0; 
  for i in range(0, n) : 
    for j in range(0, n) : 
      err = err + (B[i, j] - C[i, j])**2; 
  err; # unevaluated? 

  B = Matrix([ [ random.uniform(-1, +1)  # random matrices 
    for j in range(0, n) ] for i in range(0, n) ]); 
  C = Matrix([ [ random.uniform(-1, +1)  
    for j in range(0, n) ] for i in range(0, n) ]); 
  
  err = 0; 
  for i in range(0, n) : 
    for j in range(0, n) : 
      err = err + (B[i, j] - C[i, j])**2; 
  err; # evaluated! 
# end if 

################################################################################

