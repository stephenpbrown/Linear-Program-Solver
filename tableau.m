function tableau(A,b,c)
  
  costchoice = input("For most negative reduced cost, choose 0, for smallest index with negative reduced cost choose 1: ");

  leavechoice = input("Enter 0 for smallest index rule, enter 1 for lexicographical rule: ");
 
  tic; # Start timer after choice

  [m,n] = size(A);
  format rat compact;

  # Scale b to make sure all b(i) are positive
  newb = [];
  for i = 1:length(b)
    if(b(i) < 0) newb(i) = -1*b(i);
    else newb(i) = b(i);
    endif
  end
  b = newb';

  M = 10000;
  T = [-sum(b)*M [c' M*ones(1,m)]; b A eye(m)];
  c = [c; M*ones(m,1)];
  basis = find(c == M);
  OB = basis; # Save original basis to check for infeasibility
  
  Binv = eye(m);
  cB = c(basis);
  cp = c' - cB'*Binv*[A eye(m)];
  T(1,2:end) = cp;
  
  [m,n] = size(T);

  tol = 10^(-6);
  max_iter = 1000;

  for iter = 1:max_iter
    # Find all reduced costs
    c_r = T(1, 2:end);
    # T
    # Check if solution has been found
    J = find(c_r < -tol);
    if(isempty(J))
	  # Check if any row needs to be scaled corresponding to the ending basis
	  for i = 1:length(basis)
	    for j = 2:m
          if(T(j,basis(i)+1) != 1 && T(j,basis(i)+1) != 0) 
		    T(j,:) = T(j,:)/T(j,basis(i)+1);
          end
        end
      end

      disp("\nSolution found!\n");
      fprintf("Number of Iterations: %d\n\n", iter); 
      for i = 1:length(basis)
        fprintf("x%d = %f\n",basis(i),T(i+1,1)); 
      end
      c(c == M) = 0;
	  z = -c(basis)'*T(2:end,1);
      fprintf("\nz* = %f\n\n", z);

      # Find infeasibility
	  infeasible = false;
	  for i = 1:length(basis)
        if(!isempty(find(basis(i) == OB)))
			if(T(i+1,1) > 0)
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
   
    # Choose the input based on choosing the most negative reduced cost
    if(costchoice == 0) # a: Most negative reduced cost enters
      j = find(c_r == min(c_r));
      j = j(1)+1;
    end
    
    if(costchoice == 1) # b: Negative reduced cost with smallest index enters
      j = find(c_r < -tol);
      j = j(1)+1;
    end

	# Check for unbounded solution	
	u = T(2:end,j);
	I = find(u > 0);
	if(isempty(I) == 1)
		disp("Unbounded solution");
		toc # End timer
		return
	end

    # Lexicographical rule
    if(leavechoice == 1)
      # Set all possible pivots to 1
      inds = find(T(2:m,j)>0)+1; # inds is all the positive rows
      ratios=T(inds,1)./T(inds,j);
      theta=min(ratios);
      candxjs = basis(inds(find(ratios==theta))-1);
      xl = min(candxjs);

      for i = 1:length(inds)
        T(inds(i),:) = T(inds(i),:)/T(inds(i),j); # T(row,col)
      end

      val = [];
      for i = 1:n
        val = find(T(inds,i) == min(T(inds,i))); # Returns the indice of the minumum value
        if(length(val) == 1) # Makes sure there are no ties
          l = inds(val);
          break
        end  
        inds = inds(val); # Exclude non tied rows by updating possible indices
      end

      for i=setdiff(1:m,l)
        T(i,:) = T(i,:) - T(i,j)*T(l,:);
      end

      basis(find(basis==xl)) = j-1;
    end # End of lexicographical rule
    
    
    # Smallest index rule
    if(leavechoice == 0)    
      # Pick the smallest index
      inds = find(T(2:m,j)>0)+1; 
      ratios=T(inds,1)./T(inds,j); 
      theta=min(ratios);
      candxjs = basis(inds(find(ratios==theta))-1);

	  ## Changed these two lines ##
      xl = min(candxjs); # Finds the minimum candidate to use
      l = find(basis == xl)+1; # Finds the proper index to use for the leaving row
	  ##    Split up xl and l    ##  

      T(l,:) = T(l,:)*1/T(l,j);
      for i=setdiff(1:m,l)
        T(i,:) = T(i,:) - T(i,j)*T(l,:);
      end
    
      # Update basis
      basis(find(basis==xl)) = j-1;
    end # End of smallest index rule
    
  end # End of the main for loop
  fprintf("\nMax number of iterations reached (%d)! Cycling?\n\n",iter);
  toc # End timer
  return
end # End of the function
