#' fluspect
#'
#' \code{fluspect} calculates reflectance and transmittance spectra of a leaf using FLUSPECT-B,
#' plus four excitation-fluorescence matrices
#'
#' More information: \href{https://doi.org/10.1016/j.rse.2016.09.017}{Fluspect-B: A model for leaf fluorescence, reflectance and transmittance spectra. Vilfan et al., 2016}\cr\cr
#' Original version in MatLab: \href{https://github.com/Christiaanvandertol/Fluspect}{github.com/Christiaanvandertol/Fluspect}
#'
#' @author Nastassia Vilfan, Christiaan van der Tol, Onno Muller, Uwe Rascher, Wouter Verhoef (Original version in Matlab)
#' @author Alberto Hornero (Ported version into R)
#'
#' @param leafbio Data Frame. It contains: Cab, Cca, Cw, Cdm, Cs, N, fqe
#' @param spectral List. (Optional) A spectral object created with \link{define.bands}. A default spectral object is used
#' if the user does not indicate any.
#' @param optipar Data Frame. (Optional) It contains: nr, Kdm, Kab, Kca, Kw, Ks, phiI, phiII. A default optipar object is used
#' if the user does not indicate any.
#'
#' @return a list which contains:
#' * refl (reflectance)
#' * tran (transmittance)
#' * Mb (backward scattering fluorescence matrix, I for PSI and II for PSII)
#' * Mf (forward scattering fluorescence matrix,  I for PSI and II for PSII)
#'
#' @importFrom pracma expint
#'
#' @export
#'
#' @examples
#' leafbio <- data.frame(Cab = 70, Cca = 30, Cw = 0.013, Cdm = 0.024, Cs = 0.0, N = 4.09, fqe = 0.02)
#' leafopt <- fluspect(leafbio)
#' plot(leafopt$refl)
#'
#' @md
fluspect <- function(leafbio, spectral = define.bands(), optipar = NULL) {

  if(is.null(optipar))
    optipar = db.optipar

  leafopt <- list()

  ## Internal function which calculates TAV
  calc.tav <- function(alfa, nr) {

    rd <- pi/180
    n2 <- nr^2
    np <- n2+1
    nm <- n2-1
    a  <- (nr+1)*(nr+1)/2
    k  <- -(n2-1)*(n2-1)/4
    sa <- sin(alfa*rd)

    b1 = 0
    if(alfa!=90)
      b1 = sqrt((sa^2-np/2)*(sa^2-np/2)+k)
    b2 <- sa^2-np/2
    b  <- b1-b2
    b3 <- b^3
    a3 <- a^3
    ts <- (k^2/(6*b3)+k/b-b/2)-(k^2/(6*a3)+k/a-a/2)

    tp1 <- -2*n2*(b-a)/(np^2)
    tp2 <- -2*n2*np*log(b/a)/(nm^2)
    tp3 <- n2*(1/b-1/a)/2
    tp4 <- 16*n2^2*(n2^2+1)*log((2*np*b-nm^2)/(2*np*a-nm^2))/(np^3*nm^2)
    tp5 <- 16*n2^3*(1/(2*np*b-nm^2)-1/(2*np*a-nm^2))/(np^3)
    tp  <- tp1+tp2+tp3+tp4+tp5
    tav <- (ts+tp)/(2*sa^2)

    return(tav)
  }

  ## parameters
  ndub <- 15 # number of doublings applied

  # Fluspect parameters
  Cab <- leafbio$Cab
  Cca <- leafbio$Cca
  Cw  <- leafbio$Cw
  Cdm <- leafbio$Cdm
  Cs  <- leafbio$Cs
  N   <- leafbio$N
  fqe <- c(leafbio$fqe/5, leafbio$fqe)

  nr    <- optipar$nr
  Kab   <- optipar$Kab
  Kca   <- optipar$Kca
  Ks    <- optipar$Ks
  Kw    <- optipar$Kw
  Kdm   <- optipar$Kdm
  phiI  <- optipar$phiI
  phiII <- optipar$phiII

  ## PROSPECT calculations
  Kall <- (Cab*Kab + Cca*Kca + Cdm*Kdm + Cw*Kw  + Cs*Ks)/N # Compact leaf layer

  j          <- which(Kall > 0) # Non-conservative scattering (normal case)
  t1         <- (1-Kall)*exp(-Kall)
  t2         <- Kall^2*pracma::expint(Kall)
  tau        <- rep(1, length(t1))
  tau[j]     <- t1[j]+t2[j]
  kChlrel    <- rep(0, length(t1))
  kChlrel[j] <- (Cab*Kab)[j]/(Kall[j]*N)

  talf <- calc.tav(59, nr)
  ralf <- 1-talf
  t12  <- calc.tav(90, nr)
  r12  <- 1-t12
  t21  <- t12/(nr^2)
  r21  <- 1-t21

  # top surface side
  denom <- 1-r21*r21*tau^2
  Ta    <- talf*tau*t21/denom
  Ra    <- ralf+r21*tau*Ta

  # bottom surface side
  t <- t12*tau*t21/denom
  r <- r12+r21*tau*t

  # Stokes equations to compute properties of next N-1 layers (N real)
  # Normal case

  D  <- sqrt((1+r+t)*(1+r-t)*(1-r+t)*(1-r-t))
  rq <- r^2
  tq <- t^2
  a  <- (1+rq-tq+D)/(2*r)
  b  <- (1-rq+tq+D)/(2*t)

  bNm1   <- b^(N-1)
  bN2    <- bNm1^2
  a2     <- a^2
  denom  <- a2*bN2-1
  Rsub   <- a*(bN2-1)/denom
  Tsub   <- bNm1*(a2-1)/denom

  # Case of zero absorption
  j       <- which(r+t >= 1)
  Tsub[j] <- t[j]/(t[j]+(1-t[j])*(N-1))
  Rsub[j] <- 1-Tsub[j]

  # Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
  denom <- 1-Rsub*r
  tran  <- Ta*Tsub/denom
  refl  <- Ra+Ta*Rsub*t/denom

  leafopt$refl    <- refl
  leafopt$tran    <- tran
  leafopt$kChlrel <- kChlrel

  # From here a new path is taken: The doubling method used to calculate
  # fluoresence is now only applied to the part of the leaf where absorption
  # takes place, that is, the part exclusive of the leaf-air interfaces. The
  # reflectance (rho) and transmittance (tau) of this part of the leaf are
  # now determined by "subtracting" the interfaces

  Rb  <- (refl-ralf)/(talf*t21+(refl-ralf)*r21) # Remove the top interface
  Z   <- tran*(1-Rb*r21)/(talf*t21) # Derive Z from the transmittance

  rho <- (Rb-r21*Z^2)/(1-(r21*Z)^2) # Reflectance and transmittance
  tau <- (1-Rb*r21)/(1-(r21*Z)^2)*Z # of the leaf mesophyll layer
  t   <- tau
  r   <- rho
  r[rho<0] <- 0 # Avoid negative r

  # Derive Kubelka-Munk s and k

  I_rt     <- (r+t)<1
  D[I_rt]  <-  sqrt((1 + r[I_rt] + t[I_rt]) *
                      (1 + r[I_rt] - t[I_rt]) *
                      (1 - r[I_rt] + t[I_rt]) *
                      (1 - r[I_rt] - t[I_rt]))
  a[I_rt]  <- (1 + r[I_rt]^2 - t[I_rt]^2 + D[I_rt]) / (2*r[I_rt])
  b[I_rt]  <- (1 - r[I_rt]^2 + t[I_rt]^2 + D[I_rt]) / (2*t[I_rt])
  a[!I_rt] <- 1
  b[!I_rt] <- 1

  s        <- r/t
  I_a      <- (a>1 & a!=Inf)
  s[I_a]   <- 2*a[I_a] / (a[I_a]^2 - 1) * log(b[I_a])

  k        <- log(b)
  k[I_a]   <- (a[I_a]-1) / (a[I_a]+1) * log(b[I_a])
  kChl     <- kChlrel * k

  ## Fluorescence of the leaf mesophyll layer
  if (leafbio$fqe > 0) {
    wle <- spectral$wlE # excitation wavelengths, transpose to column
    wlf <- spectral$wlF # fluorescence wavelengths, transpose to column
    wlp <- spectral$wlP # PROSPECT wavelengths, kept as a row vector

    minwle <- min(wle)
    maxwle <- max(wle)
    minwlf <- min(wlf)
    maxwlf <- max(wlf)

    # indices of wle and wlf within wlp
    Iwle <- which(wlp>=minwle & wlp<=maxwle)
    Iwlf <- which(wlp>=minwlf & wlp<=maxwlf)

    eps <- 2^(-ndub)

    # initialisations
    te <- 1-(k[Iwle]+s[Iwle]) * eps
    tf <- 1-(k[Iwlf]+s[Iwlf]) * eps
    re <- s[Iwle] * eps
    rf <- s[Iwlf] * eps

    sigmoid <- 1/(1+ outer(exp(-wlf/10), exp(wle/10))) # matrix computed as an outer product
    # plot(sigmoid[,1], type = 'l', ylim = c(0,1), col = sample(colours(), 1))
    # for(i in 2:ncol(sigmoid))
    #   lines(sigmoid[,i], col = sample(colours(), 1) )

    # Other factor .5 deleted, since these are the complete efficiencies
    # for either PSI or PSII, not a linear combination

    # matrix multiplication preceeds element-wise multiplication
    MfI  <- MbI  <- fqe[1] * ((.5*phiI[Iwlf])*eps) %*% t(as.matrix(kChl[Iwle])) * sigmoid
    MfII <- MbII <- fqe[2] * ((.5*phiII[Iwlf])*eps) %*% t(as.matrix(kChl[Iwle])) * sigmoid

    Ih <- matrix(1, 1, length(te)) # row of ones   ## rep(1, length(te))
    Iv <- matrix(1, length(tf), 1) # column of ones

    # Doubling routine
    for (i in 1:ndub) {
      xe <- te/(1-re*re);  ten <- te*xe;  ren <- re*(1+ten)
      xf <- tf/(1-rf*rf);  tfn <- tf*xf;  rfn <- rf*(1+tfn)

      A11  <- xf %*% Ih + Iv %*% xe;              A12 <- (xf %*% t(xe)) * (rf %*% Ih + Iv %*% re)
      A21  <- 1+(xf %*% t(xe)) * (1+rf%*%t(re));  A22 <- (xf*rf)%*%Ih+Iv%*%(xe*re)

      MfnI  <- MfI  * A11 + MbI  * A12
      MbnI  <- MbI  * A21 + MfI  * A22
      MfnII <- MfII * A11 + MbII * A12
      MbnII <- MbII * A21 + MfII * A22

      te   <- ten;  re  <- ren;   tf   <- tfn;   rf   <- rfn
      MfI  <- MfnI; MbI <- MbnI;  MfII <- MfnII; MbII <- MbnII
    }

    # Here we add the leaf-air interfaces again for obtaining the final
    # leaf level fluorescences.
    g1 <- MbI; g2 <- MbII; f1 <- MfI; f2 <- MfII
    Rb <- rho + tau^2 * r21 / (1-rho * r21)

    Xe <- Iv %*% (talf[Iwle] / (1-r21[Iwle] * Rb[Iwle]))
    Xf <- t21[Iwlf]/(1-r21[Iwlf]*Rb[Iwlf]) %*% Ih
    Ye <- Iv %*% (tau[Iwle]*r21[Iwle]/(1-rho[Iwle]*r21[Iwle]))
    Yf <- tau[Iwlf]*r21[Iwlf]/(1-rho[Iwlf]*r21[Iwlf]) %*% Ih

    A <- Xe * (1 + Ye*Yf) * Xf
    B <- Xe * (Ye + Yf) * Xf

    g1n <- A * g1 + B * f1
    f1n <- A * f1 + B * g1
    g2n <- A * g2 + B * f2
    f2n <- A * f2 + B * g2

    leafopt$MbI  <- g1n
    leafopt$MbII <- g2n
    leafopt$MfI  <- f1n
    leafopt$MfII <- f2n
  }

  return(leafopt)
}

#' define.bands
#'
#' \code{define.bands} defines the spectral regions for the Fluspect-B model
#'
#' Define spectral regions for SCOPE v_1.40 \cr
#' All spectral regions are defined here as row vectors \cr
#' WV Jan. 2013
#'
#' @author Nastassia Vilfan, Christiaan van der Tol, Onno Muller, Uwe Rascher, Wouter Verhoef (Original version in Matlab)
#' @author Alberto Hornero (Ported version into R)
#'
#' @export
#'
#' @return a spectral object.
#'
#' @examples
#' spectral <- define.bands()
#'
define.bands <- function() {
    # 3 spectral regions for SCOPE
    reg1 <- as.integer(seq(400, 2400, 1))
    reg2 <- as.integer(seq(2500, 15000, 100))
    reg3 <- as.integer(seq(16000, 50000, 1000))
    wlS  <- c(reg1, reg2, reg3)

    spectral <- list()
    spectral$wlS  <- wlS

    # Other spectral (sub)regions
    spectral$wlP   <- reg1 # PROSPECT data range
    spectral$wlE   <- as.integer(seq(400, 750, 1))
    spectral$wlF   <- as.integer(seq(640, 850, 1)) # chlorophyll fluorescence in E-F matrix
    spectral$wlO   <- reg1 # optical part
    spectral$wlT   <- c(reg2, reg3) # thermal part

    spectral$wlPAR <- as.integer(seq(400, 700))

    # Data used by aggreg routine to read in MODTRAN data
    spectral$SCOPEspec$nreg  <- 3
    spectral$SCOPEspec$start <- c(400, 2500, 16000)
    spectral$SCOPEspec$end   <- c(2400, 15000, 50000)
    spectral$SCOPEspec$res   <- c(1, 100, 1000)

    return(spectral)
}

#' write.leafopt
#'
#' \code{write.leafopt} writes the leafopt object as text files
#'
#' It always writes the text files in UNIX format, under the specified output path. It will override if files already exists.
#' The output filenames are:
#' * leafoptrefl.txt
#' * leafopttran.txt
#' * leafoptkChlrel.txt
#' * leafoptMbI.txt
#' * leafoptMbII.txt
#' * leafoptMfI.txt
#' * leafoptMbII.txt
#'
#' @author Alberto Hornero
#'
#' @param leafopt List. Spectral object generated with \link{fluspect}
#' @param path String. Output path. If the directory does not exist, a new one will be created.
#' @param digits Integer. (Optional) Number of digits (by default 5)
#'
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' leafbio <- data.frame(Cab = 70, Cca = 30, Cw = 0.013, Cdm = 0.024, Cs = 0.0, N = 4.09, fqe = 0.02)
#' leafopt <- fluspect(leafbio)
#' write.leafopt(leafopt, path = file.path(tempdir(), 'output'))
#'
#' @md
write.leafopt <- function(leafopt, path, digits = 5) {

  if(!dir.exists(path)) {
    dir.create(path = file.path(path), recursive = TRUE)
  }

  conn = file(paste0(path, '/leafoptrefl.txt'), 'wb')
  write(x = signif(leafopt$refl, digits = digits), file = conn, ncolumns = 1, append = FALSE)
  close(conn)

  conn = file(paste0(path, '/leafopttran.txt'), 'wb')
  write(x = signif(leafopt$tran, digits = digits), file = conn, ncolumns = 1, append = FALSE)
  close(conn)

  conn = file(paste0(path, '/leafoptkChlrel.txt'), 'wb')
  write(x = signif(leafopt$kChlrel, digits = digits), file = conn, ncolumns = 1, append = FALSE)
  close(conn)

  conn = file(paste0(path, '/leafoptMbI.txt'), 'wb')
  write.table(x = signif(leafopt$MbI, digits = digits), file = conn, append = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  close(conn)

  conn = file(paste0(path, '/leafoptMbII.txt'), 'wb')
  write.table(x = signif(leafopt$MbII, digits = digits), file = conn, append = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  close(conn)

  conn = file(paste0(path, '/leafoptMfI.txt'), 'wb')
  write.table(x = signif(leafopt$MfI, digits = digits), file = conn, append = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  close(conn)

  conn = file(paste0(path, '/leafoptMfII.txt'), 'wb')
  write.table(x = signif(leafopt$MfII, digits = digits), file = conn, append = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  close(conn)
}
