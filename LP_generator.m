function LP_generator(n, m, p)
  
  # Debug values listed in order of occurence.
  # 1 = display values, 0 = no debug output.
  
  param =       0; # Display the value of the parameters you called LP_generator with.
  debug_x =     0; # Displays the x matrix.
  debug_0A =    0; # Displays the zeroed A matrix.
  debug_h =     0; # Displays the h matrix.
  debug_d =     0; # Displays the d matrix.
  debug_c =     0; # Displays the c matrix.
  debug_t =     0; # Displays the t matrix.
  debug_A =     0; # Displays the A matrix.
  debug_b =     0; # Displays the b matrix.
  debug_newC =  0; # Displays the newC matrix (c matrix with just the costs).
  debug_col =   0; # Displays the number of columns of A and c.
  
  # Chops of decimals.
  n = fix(n); m = fix(m); p = fix(p); 
  
  # Prints parameters.
  if(param)
    printf("\nParameters: n = %d; m = %d; p = %d\n\n", n,m,p);
  endif
  
  # Makes sure parameters are greater than 0.
  if(n < 1 || m < 1 || p < 1)
    printf("Parameters cannot be less than 1.\n\n");
    return;
  endif
  
  #################################
  # Builds array of x structures. #
  #################################
  y = 1; # indexes the array.
  for ix = 1:n
    for jx = 1:m
      for kx = 1:p
        x(y).i = ix;
        x(y).j = jx;
        x(y).k = kx;
        y += 1;
      endfor
    endfor
  endfor

  # Prints out x
  if(debug_x)
    printf("===========\n X matrix.\n===========\n");
    for y = 1:(n*m*p)
      printf("x(%d).i = %d\n",y, x(y).i);
      printf("x(%d).j = %d\n",y, x(y).j);
      printf("x(%d).k = %d\n\n",y, x(y).k);
    endfor
  endif
  
  ##########################
  # Builds zeroed matrix A #
  ##########################
  rowA = (n*p)+(m*p); 
  colA = (n*m*p);
  A = zeros(rowA,colA);
  
  # Prints out zeroed A matrix.
  if(debug_0A)
    disp("=================\n Zeroed A matrix.\n=================\n")
    disp(A);
    disp("");
  endif
  
  ################################
  # Build array of h structures. #
  ################################
  y = 1; # indexes the array.
  for jx = 1:m
    for kx = 1:p
      h(y).j = jx;
      h(y).k = kx;
      h(y).hours = fix(1000*rand(1)); # psuedo random integer generator
      y += 1;
    endfor
  endfor
  
  # Prints the h's
  if(debug_h)
    printf("============\n H matrix.\n============\n");  
    for i = 1:(m*p)
      printf("h(%d).j = %d\n", i, h(i).j);
      printf("h(%d).k = %d\n", i, h(i).k);
      printf("h(%d).hours = %d\n\n", i, h(i).hours);
    endfor
  endif
  
  ################################
  # Build array of d structures. #
  ################################ 
  y = 1;
  for ix = 1:n
    for kx = 1:p
      d(y).i = ix;
      d(y).k = kx;
      d(y).demand = fix(1000*rand(1)); # psuedo random integer generator
      y += 1;
    endfor
  endfor
  
  # Prints the d's
  if(debug_d)
    printf("===========\n D matrix.\n===========\n");
    for i = 1:(n*p)
      printf("d(%d).i = %d\n", i, d(i).i);
      printf("d(%d).k = %d\n", i, d(i).k);
      printf("d(%d).demand = %d\n\n", i, d(i).demand);
    endfor
  endif
  
  ################################ 
  # Build array of c structures. #
  ################################
  
  # Create costs based on Xijk matrix.
  for ix = 1:length(x)
    flag = 0;
    
    if(ix > 1) # Make sure there is at least one cost in c.
      # Check the c array to see if X's ik are already in c.
      for jx = 1:(ix-1)
        if(x(ix).i == c(jx).i && x(ix).k == c(jx).k) # If this index is already in the C array, copy the cost.
          c(ix).i = c(jx).i;
          c(ix).k = c(jx).k;
          c(ix).cost = c(jx).cost;
          flag = 1;
          break; # stop looking for equivalent ik's
        endif
      endfor
    endif
    
    # If a unique ik was found, generate random cost.
    if(flag == 0)
      c(ix).i = x(ix).i;
      c(ix).k = x(ix).k;
      c(ix).cost = fix(1000*rand(1));
    endif
  endfor
  
  # Prints the c's
  if(debug_c)
    printf("===========\n C matrix.\n===========\n");
    for i = 1:(n*m*p)
      printf("c(%d).i = %d\n", i, c(i).i);
      printf("c(%d).k = %d\n", i, c(i).k);
      printf("c(%d).cost = %d\n\n", i, c(i).cost);
    endfor  
  endif
  
  ################################
  # Build array of t structures. #
  ################################
  y = 1; # indexes the array.
  for ix = 1:n
    for jx = 1:m
      t(y).i = ix;
      t(y).j = jx;
      t(y).time = fix(1000*rand(1)); # psuedo random integer generator
      y += 1;
    endfor
  endfor
  
  # Prints the t's
  if(debug_t)
    printf("===========\n T matrix.\n===========\n");  
    for i = 1:(n*m)
      printf("t(%d).i = %d\n", i, t(i).i);
      printf("t(%d).j = %d\n", i, t(i).j);
      printf("t(%d).time = %d\n\n", i, t(i).time);
    endfor 
  endif
  
  ################################################
  # Fill in the t values into A. (top half of A) #\
  ################################################
  for tx = 1:length(t)
    for hx = 1:length(h)
      if(t(tx).j == h(hx).j)
        for xx = 1:length(x)
          if(t(tx).i == x(xx).i && t(tx).j == x(xx).j && x(xx).k == h(hx).k)
            A(hx,xx) = t(tx).time;
          endif
        endfor
      endif
    endfor
  endfor
  
  ######################################
  # Fill in the 1's (bottom half of A) #
  ######################################
  for dx = 1:length(d)
    for xx = 1:length(x)
      if(d(dx).i == x(xx).i && d(dx).k == x(xx).k)
        A(((m*p)+dx), xx) = 1;
      endif
    endfor
  endfor
  
  # Prints matrix A.
  if(debug_A)
    disp("=================\n Filled A matrix.\n=================\n")
    disp(A);
    disp("");
  endif
  
  ################################
  # Builds b array from h and d. #
  ################################
  y = 1; z = 1;
  for i = 1:length(h)
    b(y) = h(y).hours;
    y += 1;
  endfor
  for i = 1:length(d)
    b(y) = d(z).demand;
    y += 1; z += 1;
  endfor
  
  #############
  # Invert b. #
  #############
  b = b';
  
  # Prints the b's.
  if(debug_b)
    disp("==========\n B matrix.\n==========\n")
    disp(b);
    disp("");
  endif
  
  ###################################################
  # Create a new c array with only the cost values. #
  ###################################################
  for i = 1:length(c)
   newC(i) = c(i).cost;
  endfor

  # Prints the b's.
  if(debug_newC)
    disp("=============\n newC matrix.\n=============\n")
    disp(newC);
    disp("");
  endif
  
  ############################################
  # Display A in a copy and pastable format. #
  ############################################
  disp("");
  printf("A = [");
  for i = 1:rows(A)
    for j = 1:columns(A)
      printf(" %d", A(i,j));
    endfor
    if(i != rows(A))
      printf(";")
    else
      printf(" ];\n"); # If you want the arrays to print after you copy/paste them, remove the ; between the " "
    endif
  endfor
  
  ############################################
  # Display b in a copy and pastable format. #
  ############################################
  printf("b = [");
  for i = 1:rows(b)
    printf(" %d", b(i));
    if(i == rows(b))
      printf(" ]';\n"); # If you want the arrays to print after you copy/paste them, remove the ; between the " "
    endif
  endfor
  
  ###############################################
  # Display newC in a copy and pastable format. #
  ###############################################
  printf("c = [");
  for i = 1:columns(c)
    printf(" %d", newC(i));
    if(i == columns(c))
      printf(" ]';\n"); # If you want the arrays to print after you copy/paste them, remove the ; between the " "
    endif
  endfor
  
  ################################################
  # Display ctype in a copy and pastable format. #
  ################################################
  printf("ctype = \"");
  for i = 1:rows(A)
    printf("U");
    if(i == rows(A))
      printf("\";\n");
    endif
  endfor
  
  ##################################################
  # Display vartype in a copy and pastable format. #
  ##################################################
  printf("vartype = \"");
  for i = 1:columns(A)
    printf("C");
    if(i == columns(A))
      printf("\";\n");
    endif
  endfor
  disp("");
  
  # Make sure columns(A) == columns(c)
  if(debug_col)
    printf("\n//------ Make sure col(A) == col(c) ------//\n");
    printf("columns(A) = %d\n", columns(A));  # DEBUG
    printf("columns(c) = %d\n\n", columns(c));  # DEBUG
  endif
  
endfunction
