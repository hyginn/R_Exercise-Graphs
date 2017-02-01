# utilities.R
#
# Purpose:  Utility and convenience functions for the Graphs
#           exercises unit.
#
# Version: 1.0
#
# Date:    2016  02  25
# Author:  Boris Steipe
#
# V 1.0    First code
#
# TODO:
#

ABC.pal <- colorRampPalette(c("#9BB9EB", "#DEE6FF", "#BFCDFF",
                              "#D6BBF2", "#E56F65"), bias=0.8)


.geomSeq <- function(l, fact) {
  v <- numeric(l)
  v[1] <- 1
  for (i in 2:l) {
    v[i] <- v[i-1] * fact
  }
  return(v)
}

.syllProb <- function(p, l) {
  v <- .geomSeq(l-1, 0.93)
  v <- p * (v / sum(v))
  return(c(1-p, v))
}

ABC.rSyll <- function(n=1, pOn = 0.85, pCo = 0.85) {
  # returns N random syllables
  #
  onset <- c("", "b", "d", "f", "g", "h", "j", "k",
             "l", "m", "n", "p", "qu", "r", "s", "t",
             "v", "w", "c", "z", "bl", "br", "ch",
             "cl", "cr", "dr", "fl", "fr", "gl", "gr",
             "ph", "pl", "pr", "sh", "sk", "sl", "sm",
             "sn", "sp", "st", "sw", "th", "tr", "tw",
             "wr", "scr", "shr", "spl", "spr", "squ",
             "str", "thr")

  nuc <- c("", "a", "e", "i", "o", "u", "y", "ai",
           "au", "ay", "ea", "ee",
           "ey", "ie", "oi", "oo", "ou", "oy")

  coda <- c("", "d", "f", "g", "k", "l", "m", "n", "p",
            "r", "s", "t", "w", "x", "ck", "dd", "ds",
            "dz", "ft", "ks", "ll", "lp", "lt",
            "lx", "mb", "mm", "mp", "ms", "nf", "ng",
            "nk", "ns", "nt", "pf", "pp", "ps", "rb",
            "rd", "rf", "rg", "rm", "rn", "rp", "rr",
            "rt", "rx", "sk", "sp", "ss",
            "st", "th", "ts", "tt", "wl", "wn",
            "ws", "wt", "xt")

  syll <- character()
  for (i in 1:n) {
    s <- ""
    s <- c(s, sample(onset, 1, prob = .syllProb(pOn, length(onset))))
    s <- c(s, sample(nuc,   1, prob = .syllProb(1.0, length(nuc))))
    s <- c(s, sample(coda,  1, prob = .syllProb(pCo, length(coda))))
    syll[i] <- paste(s, collapse ="", sep ="")
  }
  return(syll)
}


# [End]
