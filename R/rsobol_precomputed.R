# Load all original functions from Owen's rsobol.R.
source("Art_Owen_Rcode.R")

# Construct once the prefix indices used to select the nested
# permutation for every dimension j and every bit level k.
#
# These indices depend only on the original Sobol bits.
# They do not depend on the random seed.
.rsobol.prefixindices <- function(thebits, m = 10, s = 5) {
  if (m < 1) {
    stop("We need m >= 1")
  }
  
  n <- 2^m
  prefix_indices <- vector("list", s)
  
  for (j in 1:s) {
    prefix_indices[[j]] <- vector("list", m)
    
    # For the first bit, all points use the root permutation.
    prefix_indices[[j]][[1]] <- rep(0, n)
    
    if (m > 1) {
      # Original Sobol bits for dimension j.
      bitmat <- thebits[, j, ]
      
      for (k in 2:m) {
        # Encode the first k - 1 original bits of every point.
        # This is the same calculation performed inside
        # Owen's .rsobol.rsobolbits().
        prefix_indices[[j]][[k]] <-
          .rsobol.bits2int(
            bitmat[, 1:(k - 1), drop = FALSE]
          )
      }
    }
  }
  
  prefix_indices
}


# Apply Owen's nested uniform scrambling to precomputed Sobol bits.
#
# This follows Owen's .rsobol.rsobolbits() function, except:
#   1. thebits is received instead of constructed;
#   2. prefix_indices is received instead of recomputed.
.rsobol.rsobolbits.precomputed <- function(
    thebits,
    prefix_indices,
    m = 10,
    s = 5,
    M = 32,
    seed = 20171215
) {
  set.seed(seed)
  
  if (m < 1) {
    stop("We need m >= 1")
  }
  
  newbits <- thebits
  n <- 2^m
  
  for (j in 1:s) {
    # Seed-dependent random binary permutations.
    # This remains inside every randomization.
    theperms <- .rsobol.getpermset2(m)
    
    for (k in 1:m) {
      indices <- prefix_indices[[j]][[k]]
      
      # Same scrambling operation as Owen's original function.
      newbits[, j, k] <-
        (thebits[, j, k] +
           theperms[[k]][1 + indices]) %% 2
    }
  }
  
  # Same treatment of bits after bit m as Owen's function.
  if (M > m) {
    newbits[, , (m + 1):M] <-
      runif(n * s * (M - m)) > 0.5
  }
  
  newbits
}


# Generate one randomized Sobol point set from the precomputed
# Sobol bits and prefix indices.
#
# This returns the same type of matrix as Owen's rsobol():
#   2^m rows and s columns.
rsobol.precomputed <- function(
    thebits,
    prefix_indices,
    m = 10,
    s = 5,
    M = 32,
    seed = 20171215
) {
  if (m > M) {
    stop(
      paste(
        sep = "",
        "Cannot deliver 2**",
        m,
        " points. Note: parameter m is log2( n )."
      )
    )
  }
  
  set.seed(seed)
  
  if (m == 0) {
    return(matrix(runif(s), nrow = 1))
  }
  
  n <- 2^m
  ans <- matrix(0, n, s)
  
  newbits <- .rsobol.rsobolbits.precomputed(
    thebits = thebits,
    prefix_indices = prefix_indices,
    m = m,
    s = s,
    M = M,
    seed = seed
  )
  
  # Same conversion to doubles as Owen's rsobol().
  for (i in 1:n) {
    for (j in 1:s) {
      ans[i, j] <-
        .rsobol.bits2unif(newbits[i, j, ])
    }
  }
  
  ans
}