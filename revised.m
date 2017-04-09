function revised(A,b,c)
  # format rat compact
  costchoice = input("Enter 0 for most reduced cost, 1 for smallest index reduced cost: ");

  tic; # Start timer after input  

  [m,n] = size(A);
  
  # Scale b to make sure all b(i) are positive
  newb = [];
  for i = 1:length(b)
    if(b(i) < 0) newb(i) = -1*b(i);
    else newb(i) = b(i);
    endif
  end
  b = newb';

  M = 10000;
  A = [A eye(m)]; # Make new A with index at the end
  c = [c; M*ones(m,1)];
  basis = find(c == M);
  OB = basis; # Save original basis to check for infeasibility
  
  # Binv = eye(m); # Need inverse to calc. reduced costs.
  cB = c(basis);
  B = A(:,basis); # Need for Binv
  Binv = inv(B);

  cp = (c' - cB'*Binv*A)'; # Start with correct reduced cost
  
  # xB creates the x values corresponding to the basis.
  xB = Binv*b; 

  # Creates a 1 by n vector full of zeros. 
  x = zeros(n,1);

  # Sets the x values based on the basis = xB
  x(basis) = xB;
  
  #
  # LOOP
  #
  tol = 10^(-6);
  max_iter = 1000;
  for iter = 1:max_iter
    
    # iter
    # calculate all the reduced costs together, i.e, the entire
    # reduced cost vector
    
    if(costchoice == 1) 
	  j = find(cp == min(cp));
	  if(!(cp(j) < -tol)) # Makes sure j didn't grab a negative zero
	    j = [];	
   	  end 
	end # Smallest index negative reduced cost
    if(costchoice == 0) j = find(cp < -tol); end # Most negative reduced cost
    # A
    # cp
	# Check for solution by seeing if there are no negative reduced costs left
    if(isempty(j))
      disp("\nSolution found!\n");
      fprintf("Number of Iterations: %d\n\n",iter);

      for i = 1:length(basis)
        fprintf("x%d = %f\n",basis(i),xB(i)); 
      end
	  cB(cB == M) = 0;
      fprintf("\nz* = %f\n\n", -cB'*xB);

      # Find infeasibility
	  infeasible = false;
	  for i = 1:length(basis)
        if(!isempty(find(basis(i) == OB)))
			if(xB(i) > 0)
				infeasible = true;
			end
		end
      end
	  if(infeasible == true)
		disp("The LP is infeasible\n");
      end

	  toc # End timer
      return 
    end

    j = j(1);
    
    # We now calculate the jth basic direction. We first find dB.
    dB = -Binv *A(:,j);

    # We now pick out the indices for which d_B is < 0.
	# Check for unbounded solution
    I_mr = find(dB < 0);
	if(isempty(I_mr) == 1)
		disp("Unbounded solution");
		toc # End timer
		return
	end
    
	ratios = -xB(I_mr)./dB(I_mr);
    thetast = min(ratios);
    l = find(ratios == min(ratios)); # Finds the lowest ratio
    l = I_mr(l);
	if(length(l) != 1) # In case of a tie, it picks the row with the smallest index
		l = min(l);
	end
    # dB
    # ratios
    # thetast
    # l
    # we initialize the new bfs as y
    y = zeros(n,1); 
    y(j)=thetast;
    y(basis) = xB + thetast*dB;

    # now update the basis by replacing basis(l) with j
    basis(l) = j; 
    
    # We call the current bfs as x (we do so at the start of each
    # iteration).
    x = y;
    xB = x(basis);
    B = A(:,basis);  
    Binv = inv(B); 
    cB = c(basis);
    cp = (c' - cB'*Binv*A)';
  
  end
  fprintf("\nMax number of iterations reached (%d)! Cycling?\n\n",iter);
  toc # End timer
  return
  
end
