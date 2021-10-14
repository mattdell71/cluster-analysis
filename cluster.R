library(cluster)
library(MASS)
library(RColorBrewer)

cols <- brewer.pal(n=4, "Set1")

# Function computing the weighted distances between elements
sq.dist <- function(d, s, sq=FALSE) {
    cols <- ncol(d)
    rows <- nrow(d)
    
    dist <- do.call(rbind, lapply(1:rows, function(i) {
        a <-  matrix(as.numeric(d[i,]), nrow=rows, ncol=cols, byrow=TRUE)
        b <-  matrix(as.numeric(s[i,]), nrow=rows, ncol=cols, byrow=TRUE)
    
        diff <- (d - a)^2
        sig <- s^2 + b^2
        dist <- (rowSums(diff/sig))  
        return(dist)
    }))

    if(sq) dist <- sqrt(dist)
    return(dist)
}

load("low.rda")

# dissimilarity matrix
diss <- as.dist(sq.dist(d.low, s.low))


# clustering: divisive and agglomerative
dn <- diana(diss)
ag <- agnes(diss, method = "ward")

# silhouette: 2,3,4 groups
sil.dn.2 <- silhouette(cutree(dn, 2), diss)
sil.dn.3 <- silhouette(cutree(dn, 3), diss)
sil.dn.4 <- silhouette(cutree(dn, 4), diss)

cl.dn.2 <- cutree(dn, 2)

sil.ag.2 <- silhouette(cutree(ag, 2), diss)
sil.ag.3 <- silhouette(cutree(ag, 3), diss)
sil.ag.4 <- silhouette(cutree(ag, 4), diss)

cl.ag.2 <- cutree(ag, 2)

# PAM  2,3,4 groups

pam2 <- pam(diss, 2)
pam3 <- pam(diss, 3)
pam4 <- pam(diss, 4)

cl.pam.2 <- pam2$clustering

# Average silhouette, all 3 methods
sm2 <- c(summary(sil.ag.2)$avg, summary(sil.dn.2)$avg, pam2$silinfo$avg)
sm3 <- c(summary(sil.ag.3)$avg, summary(sil.dn.3)$avg, pam3$silinfo$avg)
sm4 <- c(summary(sil.ag.4)$avg, summary(sil.dn.4)$avg, pam4$silinfo$avg)

SIL <- rbind(sm2,sm3,sm4)
colnames(SIL) <- c("Agnes", "Diana", "PAM")
rownames(SIL) <- 2:4
round(SIL, 2)


############### MC resampling

nrun <- 50
nd <- nrow(d.low)

res <- do.call(rbind, mclapply(1:nrun, function(i) {
    M <- do.call(rbind, lapply(1:nd, function(j) {
        SIG <- diag(as.numeric(s.low[j,]^2))
        mvrnorm(1, mu=as.numeric(d.low[j,]), Sigma=SIG)
    }))
    M <- as.matrix(M)
      
    diss <- as.dist(sq.dist(M, s.low))

    dn <- diana(diss)

    sil.dn.2 <- silhouette(cutree(dn, 2), diss)
    sil.dn.3 <- silhouette(cutree(dn, 3), diss)
    sil.dn.4 <- silhouette(cutree(dn, 4), diss)

    s2 <- summary(sil.dn.2)$avg.width
    s3 <- summary(sil.dn.3)$avg.width
    s4 <- summary(sil.dn.4)$avg.width

    return(c(s2,s3,s4))
}))

boxplot(res, outline=FALSE, names=2:4, xlab="Groups", ylab="Silhouette width")
