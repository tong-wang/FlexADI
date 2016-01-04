setwd("~/Dropbox/XCode/FlexADI.git/Homogeneous")

## read the output files
data.ADI.I <- read.table("Result_ADI_Homo_I.txt", header=TRUE)
data.ADI.II <- read.table("Result_ADI_Homo_II.txt", header=TRUE)
data.ADIF.I <- read.table("Result_ADIF_Homo_I.txt", header=TRUE)
data.ADIF.II <- read.table("Result_ADIF_Homo_II.txt", header=TRUE)

## Figure 3. Plot the optimal state-dependent (s(V), S(V)) policy
# extract the policy from output
Vs <- 0:30
policy <- data.frame(V=Vs, s=0, S=0)
for (V in Vs) {
    policy[policy$V == V, "s"] <- data.ADIF.II[1, paste0("s", V)]
    policy[policy$V == V, "S"] <- data.ADIF.II[1, paste0("S", V)]
}
policy

# generate (Y, U, V) data for plotting
Us <- -100:50
data.fig3 <- expand.grid(U=Us, V=Vs)
data.fig3$Y <- 0
for (V in Vs) {
    for (U in Us) {
        data.fig3[data.fig3$U == U & data.fig3$V == V, "Y"] <- ifelse(U <= data.ADIF.II[1, paste0("s", V)], data.ADIF.II[1, paste0("S", V)], U)
    }
}

# 3D-plot on (Y, U, V) 
library("lattice")
wireframe(Y ~ U * V, data = data.fig3, xlab="U", ylab = "V", zlab="Y", main = "Tau2_int", screen = list(z = -20, x = -75))


## Figure 4. Cost comparison between ADI and ADIF
ADI.fig4 <- rbind(data.ADI.I[, 1:8], data.ADI.II[, 1:8])
ADIF.fig4 <- rbind(data.ADIF.I[, 1:8], data.ADIF.II[, 1:8])

wireframe(Cost ~ L * T, data = ADI.fig4, xlab="L", ylab = "T", zlab="Cost", main = "Cost", screen = list(z = -20, x = -75))
wireframe(Cost ~ L * T, data = ADIF.fig4, xlab="L", ylab = "T", zlab="Cost", main = "Cost", screen = list(z = -20, x = -75))



myGrid <- cbind(ADI.fig4, Cost2=ADIF.fig4$Cost)
mypanel <- function(x,y,z,z2,...) {
    panel.wireframe(x,y,z,...)
    panel.wireframe(x,y,z2,...)
}
wireframe(Cost2 ~ - L * T, data=myGrid, xlab="L", ylab="T", zlab="Cost", zlim=c(600, 1500),
          panel=mypanel, z2=myGrid$Cost, screen = list(z = -50, x = -60), scales = list(arrows=FALSE, cex= .8, col = "black", font = 3, x = list(labels = seq(4, 0, -1)), y = list(at=0:2, labels = seq(0, 2, 1))))

