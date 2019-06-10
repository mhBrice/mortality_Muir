mozaicplot <- function(tab, col) {
  
  di <- dim(tab)
  
  prop <- function(x) x/sum(x)
  propc <- function(x) cumsum(x/sum(x))
  
  coords <- rev(expand.grid(rev(dimnames(tab))))
  names(coords) <- paste0(c("row", "col", "sli"), "names")
  
  ## relative height for each column
  hc <- apply(tab, 2, apply, 1, sum) %>%
    apply(1, prop) %>%
    apply(2, function(x) 1 - cumsum(x))
  
  coords$ybottom <- rep(as.vector(hc), each = di[3])
  coords$ytop <- rep(as.vector(rbind(1, hc[-di[2],])), each = di[3])
  
  ## width columns
  wc <- apply(tab, 1, sum) %>% propc
  pr <- prod(di[2:3])
  ba <- rep(c(0, wc[-length(wc)]), each = pr)
  le <- rep(diff(c(0,wc)), each = pr)
  v1 <- apply(tab, c(2,1), function(x) head(c(0, propc(x)), -1)) %>% as.vector
  v2 <- apply(tab, c(2,1), function(x) propc(x)) %>% as.vector
  
  coords$xleft <- ba + le*v1
  coords$xright <- ba + le*v2
  
  coords$colors <- col[as.numeric(coords$colnames)]
  
  
  plot0(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  for (i in 1:nrow(coords)) {
    rect(
      coords$xleft[i], coords$ybottom[i], coords$xright[i], coords$ytop[i],
      col = coords$colors[i], 
      density = c(30, NULL)[i%%2 + 1],
      border = "grey15", lwd = 1.5, xpd = NA
    )
  }
  abline(v = coords$xleft[c(9,17,25,33)], lwd = 2.5, col = "grey15")
  axis(2, at = c(.1, .5, .86, .98), labels = rev(unique(coords$colnames)), 
       las =1, line = -.8, tick = F, cex.axis = .75)
  axis(3, at = c(.2, .55, .72, .85, 1), labels = unique(coords$rownames), cex.axis = .75, font = 3, tick = F,
       line = -.8, xpd=NA)
  legend(1,.5, legend = c("Alive", "dead"), fill = col[4], density = c(NA,30), cex = .8,
         xpd = NA, bty = "n", border = "grey15")
  
}
