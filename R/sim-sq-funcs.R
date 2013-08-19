# sim-sq-funcs
# Functions for placing random or uniform reserve networks on a single square
# species range.  The individual reserves are also square.

# The initial set of functions are abstracted from the main program 
# Sim1Square2.R, written in 2011. The last two functions are the workhorses,
# found in Sim1SquareFuncs.R, also written in 2011. They were painstakingly
# coded, and are painful to read. I did some testing on them at the time, so
# I am trusting that they are correct.



#' Calculate the proportion of a square species range that is covered by a
#' random reserve network
#' 
#' @param res.num logical square matrix, with row/column dimension `PR.x`,
#'    indicating which of the planning units within the planning region are
#'    reserves
#' @param range.x length of one side of the species range
#' @param PR.x number of planning units along the side of the planning region
#' @param c1 length-2 vector giving the lower left coordinate of the species 
#'    range. Defaults to a random location in the planning region
#' @param c2 length-2 vector giving the upper right coordinate of the species
#'    range. Defaults to a value that will produce a square range with side
#'    range.x. \textbf{This should never be set to a non-default value!}
#' @return fraction of the species range that is included in reserves
#' @details There appears to be no error checking to ensure that the whole 
#'    species range is within the planning region. In fact, it seems that the 
#'    upper bound on c1 should be `PR.x - range.x`. In addition, the ability
#'    to specify `c2` suggests that I intended this to be used with rectangular
#'    species ranges; but the code explicitly assumes that the range has
#'    length `range.x` in both dimensions.
coverage.rand = function(res.num, range.x, PR.x, c1 = runif(2, 1, PR.x),
                         c2 = c1 + range.x) {
  cf1 = floor(c1)
  cf2 = floor(c2)
  if (range.x < 1) {
    if (cf1[1] == cf2[1]) {
      if (cf1[2] == cf2[2]) { # Range entirely w/in one PU
        frac = res.num[cf1[1],cf1[2]]
      } else { # Range in 2 PUs in direction 2
        div = (cf2[2] - c1[2])/range.x
        frac = sum(res.num[cf1[1],cf1[2]:cf2[2]] * c(div, 1-div))
      }
    } else { # spans a PU boundary in direction 1
      if (cf1[2] == cf2[2]) { # Range in 2 PUs in direction 1
        div = (cf2[1] - c1[1])/range.x
        frac = sum(res.num[cf1[1]:cf2[1],cf1[2]] * c(div, 1-div))
      } else { # Range in 4 PUs
        div1 = (cf2[1] - c1[1])/range.x
        div2 = (cf2[2] - c1[2])/range.x
        frac = sum(res.num[cf1[1]:cf2[1],cf1[2]:cf2[2]] *
                     matrix(c(div1*div2, (1-div1)*div2, div1*(1-div2), (1-div1)*(1-div2)),2,2))
      }
    }
  } else if (range.x == 1) { # Range in 4 PUs
    div1 = (cf2[1] - c1[1])/range.x
    div2 = (cf2[2] - c1[2])/range.x
    frac = sum(res.num[cf1[1]:cf2[1],cf1[2]:cf2[2]] *
                 matrix(c(div1*div2, (1-div1)*div2, div1*(1-div2), (1-div1)*(1-div2)),2,2))
  } else {
    if (cf2[1]-cf1[1] == 1 & cf2[2]-cf1[2] == 1 ) { # Range in 4 PUs
      div1 = (cf2[1] - c1[1])/range.x
      div2 = (cf2[2] - c1[2])/range.x
      frac = sum(res.num[cf1[1]:cf2[1],cf1[2]:cf2[2]] *
                   matrix(c(div1*div2, (1-div1)*div2, div1*(1-div2), (1-div1)*(1-div2)),2,2))
    } else {  # General purpose calculation
      # First do corners
      d11 = cf1[1] + 1 - c1[1]
      d12 = cf1[2] + 1 - c1[2]
      d21 = c2[1] - cf2[1]
      d22 = c2[2] - cf2[2]
      Ares = res.num[cf1[1],cf1[2]] * d11*d12 +
        res.num[cf2[1],cf1[2]] * d21*d12 +
        res.num[cf2[1],cf2[2]] * d21*d22 +
        res.num[cf1[1],cf2[2]] * d11*d22
      
      # Now do sides as needed
      if (cf2[1] - cf1[1] > 1) {
        for (ii in (cf1[1]+1):(cf2[1]-1)) {
          Ares = Ares + res.num[ii,cf1[2]] * d12 + res.num[ii,cf2[2]] * d22
        }
      }
      if (cf2[2] - cf1[2] > 1) {
        for (ii in (cf1[2]+1):(cf2[2]-1)) {
          Ares = Ares + res.num[cf1[1],ii] * d11 + res.num[cf1[2],ii] * d12
        }
      }
      # Now do interior PUs as needed
      if (cf2[1] - cf1[1] > 1 & cf2[2] - cf1[2] > 1) {
        Ares = Ares + sum(res.num[(cf1[1]+1):(cf2[1]-1), (cf1[2]+1):(cf2[2]-1)])
      }
      frac = Ares/range.x^2
    }
  }
  return(frac)
}

coverage.unif = function(x.unif, rep.unif, n.unif, range.x,
                         c1 = sort(runif(2, 0, rep.unif)), c2 = c1 + range.x) {
  ijmax = ceiling(c2/rep.unif)
  Ares = 0
  for (ii in 1:ijmax[1]) {
    ibase = (ii-1) * rep.unif
    for (jj in 1:ijmax[2]) {
      jbase = (jj-1) * rep.unif
      Ares = Ares + max(0, min(ibase+x.unif, c2[1]) - max(ibase, c1[1])) *
        max(0, min(jbase+x.unif, c2[2]) - max(jbase, c1[2]))
    }
  }
  prop = Ares/range.x^2
  return(prop)
}