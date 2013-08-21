#' Calculate the proportion of a square species range that is covered by a
#' uniform reserve network
#' 
#' @param x.unif the length of one side a of a reserve
#' @param rep.unif the linear dimension of a reserve+gap unit
#' @param n.unif unused parameter(!)
#' @param range.x length of one side of the species range
#' @param c1 length-2 vector giving the lower left coordinate of the species 
#'    range. Defaults to a random location in the planning region
#' @return fraction of the species range that is included in reserves
coverage.unif <-
function(x.unif, rep.unif, n.unif, range.x,
                         c1 = sort(runif(2, 0, rep.unif))) {
  c2 = c1 + range.x
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
