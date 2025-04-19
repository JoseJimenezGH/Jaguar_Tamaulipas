logit <- function(x) {log(x/(1-x))}

spatial.plot<- function (x, y, add = FALSE, cx = 1, col = "gray")
{
    nc <- as.numeric(cut(y, 10))
    if (!add)
        plot(x, pch = " ", asp = 1)
    if (col == "gray") {
        cc <- seq(3, 17, , 10)/20
        cc <- gray(cc)
    }
    else cc <- terrain.colors(10)
    points(x, pch = 20, col = cc[nc], cex = cx)
    image.scale(y, col = cc)
}

#Hay 5 tipos de spiderplot:
#1. spiderplot: clasico
#2. spiderplotJJ: posibilidad de buffer alrededor de las trampas
#3. spiderplotJJ2: diferentes colores por animal, diferentes grosores en los segmentos y buffer.
#4. spiderplotJJ3: spiderplot a un plot existente, con diferentes colores, grosores y buffer
#5. spiderplotJJ4: spiderplot clasico a un plot existente

spiderplot<-function (y, traplocs)
{
    dither <- FALSE
    dx <- max(traplocs[, 1]) - min(traplocs[, 1])
    dy <- max(traplocs[, 2]) - min(traplocs[, 2])
    dx <- 0.01 * dx
    dy <- 0.01 * dy
    if (length(dim(y)) == 3) {
        if (dim(y)[2] == nrow(traplocs)) {
            nind <- dim(y)[1]
            ntraps <- dim(y)[2]
            nocc <- dim(y)[3]
            newy <- array(NA, dim = c(nind, nocc, ntraps))
            for (i in 1:nind) {
                newy[i, 1:nocc, 1:ntraps] <- t(y[i, , ])
            }
            y <- newy
        }
        y3d <- y
        J <- dim(y3d)[3]
        T <- dim(y3d)[2]
        nind <- dim(y3d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y3d[i, t, ]
                if (sum(aa) > 0) {
                  aa <- traplocs[aa > 0, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            delta <- c(runif(1, -dx, dx), runif(1, -dy, dy)) *
                ifelse(dither, 1, 0)
            points(avg.s[i, 1] + delta, avg.s[i, 2] + delta,
                pch = "S", cex = 1, col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    if (length(dim(y)) == 2) {
        y2d <- y
        J <- nrow(traplocs)
        T <- dim(y2d)[2]
        nind <- dim(y2d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y2d[i, t]
                if (aa <= J) {
                  aa <- traplocs[aa, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            points(avg.s[i, 1], avg.s[i, 2], pch = "S", cex = 1,
                col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    points(traplocs, pch = 20)
    Cx <- mean(traplocs[, 1])
    Cy <- mean(traplocs[, 2])
    xcent <- sqrt((avg.s[, 1] - Cx)^2 + (avg.s[, 2] - Cy)^2)
    list(xcent = xcent, avg.s = avg.s, center = c(Cx, Cy))
}


spiderplotJJ<-function (y, traplocs, buffer)
{
    dither <- FALSE
    dx <- max(traplocs[, 1]) - min(traplocs[, 1])
    dy <- max(traplocs[, 2]) - min(traplocs[, 2])
    dx <- 0.01 * dx
    dy <- 0.01 * dy

    Xl<-min(traplocs[,1])-buffer
    Xu<-max(traplocs[,1])+buffer
    Yl<-min(traplocs[,2])-buffer
    Yu<-max(traplocs[,2])+buffer
    xlim<-c(Xl,Xu)
    ylim<-c(Yl,Yu)

    if (length(dim(y)) == 3) {
        if (dim(y)[2] == nrow(traplocs)) {
            nind <- dim(y)[1]
            ntraps <- dim(y)[2]
            nocc <- dim(y)[3]
            newy <- array(NA, dim = c(nind, nocc, ntraps))
            for (i in 1:nind) {
                newy[i, 1:nocc, 1:ntraps] <- t(y[i, , ])
            }
            y <- newy
        }
        y3d <- y
        J <- dim(y3d)[3]
        T <- dim(y3d)[2]
        nind <- dim(y3d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5, xlim=xlim, ylim=ylim)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y3d[i, t, ]
                if (sum(aa) > 0) {
                  aa <- traplocs[aa > 0, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            delta <- c(runif(1, -dx, dx), runif(1, -dy, dy)) *
                ifelse(dither, 1, 0)
            points(avg.s[i, 1] + delta, avg.s[i, 2] + delta,
                pch = "S", cex = 1, col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    if (length(dim(y)) == 2) {
        y2d <- y
        J <- nrow(traplocs)
        T <- dim(y2d)[2]
        nind <- dim(y2d)[1]
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
        avg.s <- matrix(NA, nrow = nind, ncol = 2)
        for (i in 1:nind) {
            tmp <- NULL
            for (t in 1:T) {
                aa <- y2d[i, t]
                if (aa <= J) {
                  aa <- traplocs[aa, ]
                  tmp <- rbind(tmp, aa)
                }
            }
            avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
            points(avg.s[i, 1], avg.s[i, 2], pch = "S", cex = 1,
                col = "red")
            for (m in 1:nrow(tmp)) {
                if (nrow(tmp) > 1)
                  lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i,
                    2], tmp[m, 2]))
            }
        }
    }
    points(traplocs, pch = 20)
    Cx <- mean(traplocs[, 1])
    Cy <- mean(traplocs[, 2])
    xcent <- sqrt((avg.s[, 1] - Cx)^2 + (avg.s[, 2] - Cy)^2)
    list(xcent = xcent, avg.s = avg.s, center = c(Cx, Cy))
}


spiderplotJJ2<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
            las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = clr[i], lwd = lwd)
        }
        points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
            cex = 1)
        points(mu.x[[s]], mu.y[[s]], pch = 16, cex = 1.5, col = clr)
    }
    par(op)
}

spiderplotJJ3<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        #plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
        #    las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = clr[i], lwd = lwd)
        }
        #points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
        #    cex = 1)
        points(mu.x[[s]], mu.y[[s]], pch = 16, cex = 1.5, col = clr)
    }
    par(op)
}

spiderplotJJ4<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        #plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
        #    las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = 1, lwd = lwd)
        }
        #points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
        #    cex = 1)
        points(mu.x[[s]], mu.y[[s]], pch = 16, cex = 1.5, col = 'violet')
    }
    par(op)
}



e2dist <- function (x, y) {  # Function from scrbook package to calculate the distance between locations in 2 matrices.
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


SCRdensity<-function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col = terrain.colors(ncolors))
    image.scale(Dn, col = terrain.colors(ncolors))
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}


# Con este comando sacamos un distribución real del espacio entre el eje de
# las X y el de las Y. Para que sea regular el pixel, lo ajusto con asratio
# plot(X, asp=1)
# asratio <- (Yu-Yl)/(Xu-Xl)
# nx<-75
# ny<-75* asratio

SCRdensityJJ1<-function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10, asp=1)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col=gray.colors(ncolors, start=1, end=0),asp=1)
    image.scale(Dn, col=gray.colors(ncolors, start=1, end=0))
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}




GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- lam0*exp(-d2/(2*sigma^2))
     return(ans)
    }
  }
)

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lambda = double(1),
  log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dpois(x[j], lambda[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rPoissonVector  <- nimbleFunction(
  run = function(n = integer(), lambda = double(1)) {
    J <- length(lambda)
    ans<- numeric(J)
    for(j in 1:J)
      ans[j] <- rpois(1, lambda[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dPoissonVector = list(
    BUGSdist = "dPoissonVector(lambda)",
    Rdist = "dPoissonVector(lambda)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = double(1)', 'lambda = double(1)'))
))

## sampler to jointly update y.un[1:M,j] so that they sum to n[j]
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    nnidd<-control$nnidd
    j<-control$j
    M<-control$M
    calcNodes <- model$getDependencies(target)
  },

run = function() {
  lam.curr <- model$lam[1:M,j] # individual by trap expected counts

  #Sample y[1:M,j] by reassigning n[j] using full conditional
  switch.probs <- lam.curr[1:M]/sum(lam.curr[1:M])

  #propose new ID's for nnid[j,k]
  y.latent.curr <- model$y.full[1:M,j]- model$y.obs[1:M,j]
  y.latent.prop <- rmulti(1, nnidd, switch.probs[1:M])
  model$y.full[1:M,j] <<-  model$y.obs[1:M,j] + y.latent.prop

  # initial model logProb
  model_lp_initial <- model$getLogProb(calcNodes)

  # proposal model logProb
  model_lp_proposed <- model$calculate(calcNodes)

  # log-Metropolis-Hastings ratio
  log_MH_ratio <- (model_lp_proposed + dmulti(y.latent.curr, nnidd, switch.probs, log=TRUE)) -
                  (model_lp_initial + dmulti(y.latent.prop,  nnidd, switch.probs, log=TRUE))

  # Metropolis-Hastings step
  accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)
