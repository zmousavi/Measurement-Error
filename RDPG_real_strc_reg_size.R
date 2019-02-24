
#Idea: 1 RDPG graph with latent positions (X) on H-W curve. 
#Simulate covariance(X_hat-X) entries per Avanti's theorem to generate X_hat.
#Note: we do not generate RDPG graph here.
#compute an estimate for covariance(X_hat-X)

#convert X_hat to theta_hat 
#set up graph regression with theta_hat
#Using Avanti's theorem construct an estimate for covar(X_hat-X). Find var(theta_hat)
#use the estimate to adjust measurement error 


rm(list=ls())

library("rgl")
library("MASS")
library("MCMCpack")

b0_adj_mean = c() 
b1_adj_mean =  c() 
b0_naive_mean = c() 
b1_naive_mean = c() 
b0_true_mean = c()  
b1_true_mean = c() 

b0_adj_mse = c() 
b1_adj_mse =  c() 
b0_naive_mse = c() 
b1_naive_mse = c() 
b0_true_mse = c()  
b1_true_mse = c() 
b0_mse_ratio = c()
b1_mse_ratio = c()

b0_adj_sd = c()
b1_adj_sd  = c()
b0_naive_sd = c()
b1_naive_sd = c()
b0_true_sd = c()
b1_true_sd = c()

mean_b0_adj_bias_abs = c()
mean_b1_adj_bias_abs = c()
mean_b0_naive_bias_abs = c()
mean_b1_naive_bias_abs = c()
mean_b0_true_bias_abs = c()
mean_b1_true_bias_abs = c()

g = 0
graph_size_list = c(50, 100, 500, 1000, 5000) #, 8000, 10000) 


mean_abs_X_bias_graphs =c()

for (graph_size in graph_size_list) { 
  
  g = g+1
  print(g)
  print(graph_size_list[g])
 

beta_naive = list()
beta_true = list()
beta_adj = list()
beta_true_se = list()
beta_naive_se = list()
beta_adj_se=list()
mean_abs_X_bias =c()

mc_runs = 10
r=1
for (r in c(1:mc_runs)){
  if (r%%1 ==0){
    print (r)
  }
  
  
  set.seed(r+20)
  
  #theta_true_unordered = c(0, 0.5, 1, runif(n=(graph_size-3), min=0, max=1 ))
  theta_true = seq(0, 1, by=1*(1/(graph_size-1)))
  
  theta_hw = function(th) c(th^2,2*th*(1-th),(1-th)^2)
  #theta_cube = function(th) c((th)^3,(th)^2,(th))
  X = matrix(0,nrow=graph_size,ncol=3)
  for(i in 1:graph_size) X[i,] = theta_hw((theta_true[i]))
  P = X%*%t(X)
  
  A = matrix(0,nrow=graph_size,ncol=graph_size)
  for(i in 1:(graph_size-1)) for(j in (i+1):graph_size) A[i,j] = rbinom(1,1,P[i,j])
  A = A + t(A)
  # compare: diagaug
  svdA = svd(A)
  Xhat = svdA$u %*% diag(sqrt(svdA$d))[,1:3]
  
  W_Xhat = procrustes(Xhat, X)$X.new
  
  # df=data.frame(rbind(X, Xhat,  W_Xhat))
  # df$type = c(rep(1, (graph_size)), rep(2, (graph_size)), rep(3, (graph_size)))
  # with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
  # triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")
  # 
  
  
  set.seed(r)
  beta0 = 1
  beta1 = 1
#y = b0 + b1*theta_true + rnorm(n=graph_size, 0, 1e-5)
y = beta0 + beta1*X[,1]  + rnorm(n=graph_size, 0, 1e-1)

X_true = X[,1]  #(X for regression ... not to be cofused with latent possition)
X_naive = W_Xhat[,1] 
Xz_big = cbind(rep(1, graph_size), X_true)
X_naive_big = cbind(rep(1, graph_size), X_naive)
beta_true[[r]] = solve(t(Xz_big)%*%Xz_big)%*%t(Xz_big)%*% y
beta_naive[[r]] = solve(t(X_naive_big)%*%(X_naive_big))%*%t(X_naive_big)%*% y


mean_abs_X_bias[r] = mean(abs(W_Xhat-X))

#Simulate covaraince(X_hat-X) per Avanti
Delta_inv = matrix(c(9, -9, 3, -9, 21, -9, 3, -9, 9), nrow=3)

X_noise = matrix(, ncol=3, nrow=graph_size)
#sigma_u_theta_true = c()
sigma_X1_true = c()
for (i in 1:graph_size){
  f = theta_true[i]
  S11 = (1)*((1/126) + (f/15)+(f^2)/42-(f^3)/105-2*(f^4)/35)
  S12 = (1)*((13/1260) + (4*f/105)-((f^2)/30)+(2*(f^3)/105)-((f^4)/70))
  S13 = (1)*((1/180)+(f/105)-((f^2)/70)+(f^3)/105-(f^4)/210)
  S22 = (1)*((1/45)+(4*f/105)-(2*f^2/35)+(4*f^3/105)-(2*f^4/105))
  S23 = (1)*((5/252)+(f/35)-13*(f^2)/210+(4*f^3)/105-(f^4/70))
  S33 = (1)*((2/63)+(f/7)-(73*(f^2)/210)+(5*(f^3)/21) - (2*(f^4)/35))
  S_matrix = matrix(c(S11, S12, S13, S12, S22, S23, S13, S23, S33), nrow=3, byrow=T)
  if(any(diag(S_matrix)<0)){print("yes")}
  Sigma_X = (1/graph_size)*Delta_inv%*%S_matrix%*%Delta_inv
  sigma_X1_true[i] = Sigma_X[1, 1]
  #sigma_u_theta_true[i] = (1/4)*(1*Sigma_X[1,1]-2*Sigma_X[1,3]+Sigma_X[3,3])
  #only if simulating: X_noise[i,1:3] = X[i,1:3]  + mvrnorm(n=1, rep(0, 3), Sigma=Sigma_X)
}



Z = cbind(y, X_naive_big)

#Sigma_U_theta=sigma_u_theta_est_hw
Sigma_U=sigma_X1_true

M = matrix(0,nrow=3,ncol=3)
for (i in (1:graph_size)){
  M = (Z[i,]%o%Z[i,]-diag(c(0, 0, Sigma_U[i]))) + M
}
M_z = (1/graph_size)*M
beta_adj[[r]] = solve(M_z[2:3,2:3])%*%M_z[2:3,1]

beta_true_se[[r]] = (beta_true[[r]]-c(beta0, beta1))^2
beta_adj_se[[r]] = (beta_adj[[r]]-c(beta0, beta1))^2
beta_naive_se[[r]] = (beta_naive[[r]]-c(beta0, beta1))^2

}

b0_true = c()
b1_true = c()
b0_naive = c()
b1_naive = c()
b0_adj = c()
b1_adj = c()
b0_true_se = c()
b1_true_se = c()
b0_naive_se = c()
b1_naive_se = c()
b0_adj_se = c()
b1_adj_se = c()


for (z in 1:mc_runs){
  
  
  
  b0_true[z] =  beta_true[[z]][1]
  b1_true[z] =  beta_true[[z]][2]
  
  b0_true_se[z] =  beta_true_se[[z]][1]
  b1_true_se[z] =  beta_true_se[[z]][2]
  
  b0_naive[z] = beta_naive[[z]][1]
  b1_naive[z] = beta_naive[[z]][2]
  
  b0_naive_se[z] = beta_naive_se[[z]][1]
  b1_naive_se[z] = beta_naive_se[[z]][2]
  
  b0_adj_se[z] = beta_adj_se[[z]][1]
  b1_adj_se[z] = beta_adj_se[[z]][2]
  
  b0_adj[z] = beta_adj[[z]][1]
  b1_adj[z] = beta_adj[[z]][2]
  
  
}


b0_adj_mean[g] = mean(b0_adj) 
b1_adj_mean[g] = mean(b1_adj) 
b0_naive_mean[g] = mean(b0_naive) 
b1_naive_mean[g] = mean(b1_naive) 
b0_true_mean[g] = mean(b0_true) 
b1_true_mean[g] = mean(b1_true)
mean_b0_adj_bias_abs[g] = mean(abs(b0_adj-beta0))
mean_b1_adj_bias_abs[g] =  mean(abs(b1_adj-beta1))
mean_b0_naive_bias_abs[g] = mean(abs(b0_naive-beta0))
mean_b1_naive_bias_abs[g] = mean(abs(b1_naive-beta1))
mean_b0_true_bias_abs[g] = mean(abs(b0_true-beta0))
mean_b1_true_bias_abs[g] = mean(abs(b1_true-beta1))

b0_adj_sd[g] = sd(b0_adj) 
b1_adj_sd[g] = sd(b1_adj) 
b0_naive_sd[g] = sd(b0_naive) 
b1_naive_sd[g] = sd(b1_naive) 
b0_true_sd[g] = sd(b0_true) 
b1_true_sd[g] = sd(b1_true)

b0_adj_mse[g] = mean(b0_adj_se) 
b1_adj_mse[g] = mean(b1_adj_se) 
b0_naive_mse[g] = mean(b0_naive_se) 
b1_naive_mse[g] = mean(b1_naive_se) 
b0_true_mse[g] = mean(b0_true_se) 
b1_true_mse[g] = mean(b1_true_se)

b0_mse_ratio = b0_adj_mse/b0_naive_mse 
b1_mse_ratio = b1_adj_mse/b1_naive_mse

mean_abs_X_bias_graphs[g] = mean(mean_abs_X_bias)




h = b1_true - beta1
(mean(h))^2+(sd(h))^2
b1_true_mse

f = b1_naive - beta1
(mean(f))^2+(sd(f))^2
b1_naive_mse
k = b1_adj - beta1
(mean(k))^2+(sd(k))^2
b1_adj_mse
((mean(k))^2+(sd(k))^2)/((mean(f))^2+(sd(f))^2)

setwd("~/Desktop/Feb22")

file0 = paste("beta0_real_strc_", graph_size_list[g], ".pdf", sep="")
file1 = paste("beta1_real_strc_", graph_size_list[g], ".pdf", sep="")



dev.new()
pdf(file0, 7, 5)
boxplot(b0_true, b0_naive,  b0_adj,  notch=TRUE,
        main=bquote(paste("b0 Estimate","\n graph_size:", graph_size , "\n mc_runs:",length(b1_true))),
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = beta0, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf(file1, 7, 5)
boxplot(b1_true, b1_naive,  b1_adj,  notch=TRUE,
        main=bquote(paste("b1 Estimate","\n graph_size:", graph_size , "\n mc_runs:",length(b1_true))),
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = beta1, col="red", lty=2, lwd=0.5)
dev.off()


}

min1 = min(b1_true_mean, b1_adj_mean, b1_naive_mean)
max1 = max(b1_true_mean, b1_adj_mean, b1_naive_mean)

dev.new()
pdf("beta1_size.pdf", 4, 4)
plot(graph_size_list, b1_true_mean, col="blue", ylim=c(min1, max1), main = "", xlab="", ylab="")
lines(graph_size_list, b1_true_mean, col="blue")
grid()
par(new=T)
plot(graph_size_list, b1_adj_mean, col="green", ylim=c(min1, max1), main = "", xlab="", ylab="")
lines(graph_size_list, b1_adj_mean, col="green")
par(new=T)
plot(graph_size_list, b1_naive_mean, col="red", ylim=c(min1, max1), main = "", xlab="", ylab="")
lines(graph_size_list, b1_naive_mean, col="red")
dev.off()


min0 = min(b0_true_mean, b0_adj_mean, b0_naive_mean)
max0 = max(b0_true_mean, b0_adj_mean, b0_naive_mean)
dev.new()
pdf("beta0_size.pdf", 4, 4)
plot(graph_size_list, b0_true_mean, col="blue", ylim=c(min0, max0), main = "", xlab="", ylab="")
lines(graph_size_list, b0_true_mean, col="blue")
grid()
par(new=T)
plot(graph_size_list, b0_adj_mean, col="green", ylim=c(min0, max0), main = "", xlab="", ylab="")
lines(graph_size_list, b0_adj_mean, col="green")
par(new=T)
plot(graph_size_list, b0_naive_mean, col="red", ylim=c(min0, max0), main = "", xlab="", ylab="")
lines(graph_size_list, b0_naive_mean, col="red")
dev.off()

min2 = min(b1_mse_ratio, b0_mse_ratio)
max2 = max(b1_mse_ratio, b0_mse_ratio)
dev.new()
pdf("beta_real_strc_mse_size.pdf", 4, 4)
plot(graph_size_list, b1_mse_ratio, col="blue", ylim=c(min2, max2), main = "MSE Ratio", xlab="", ylab="")
lines(graph_size_list, b1_mse_ratio, col="blue")
grid()
par(new=T)
plot(graph_size_list, b0_mse_ratio, col="green", ylim=c(min2, max2), main = "", xlab="", ylab="")
lines(graph_size_list, b0_mse_ratio, col="green")
legend("topright", legend=c("b0", "b1"), col=c("green", "blue"), lty=1, cex=0.8)
dev.off()


bias_b1_adj_2 = (b1_adj_mean-beta1)^2
bias_b1_naive_2 = (b1_naive_mean-beta1)^2
min5 = min(bias_b1_adj_2, bias_b1_naive_2)
max5 = max(bias_b1_adj_2, bias_b1_naive_2)
dev.new()
pdf("beta1_bias2_real_strc_size.pdf", 4, 4)
plot(graph_size_list, bias_b1_adj_2, col="green", ylim=c(min5, max5), main = "Bias_Squared", xlab="", ylab="")
lines(graph_size_list, bias_b1_adj_2, col="green")
grid()
par(new=T)
plot(graph_size_list, bias_b1_naive_2, col="red", ylim=c(min5, max5), main = "", xlab="", ylab="")
lines(graph_size_list, bias_b1_naive_2, col="red")
legend("topright", legend=c("b1_adj", "b1_naive"), col=c("green", "red"), lty=1, cex=0.8)
dev.off()

bias_b0_adj = (b0_adj_mean-beta0)
bias_b0_naive = (b0_naive_mean-beta0)
bias_b1_adj = (b1_adj_mean-beta1)
bias_b1_naive = (b1_naive_mean-beta1)
delta_bias_b0 = (bias_b0_naive-bias_b0_adj)/bias_b0_naive
delta_bias_b1 = (bias_b1_naive-bias_b1_adj)/bias_b1_naive
min6 = min(delta_bias_b0, delta_bias_b1)
max6 = max(delta_bias_b0, delta_bias_b1)
dev.new()
pdf("delta_bias_real_strc_size.pdf", 4, 4)
plot(graph_size_list, delta_bias_b0, col="green", ylim=c(min6, max6), main = "(Bias(Naive)-Bias(Adj))/Bias(Naive)", xlab="", ylab="")
lines(graph_size_list, delta_bias_b0, col="green")
grid()
par(new=T)
plot(graph_size_list, delta_bias_b1, col="blue", ylim=c(min6, max6), main = "", xlab="", ylab="")
lines(graph_size_list, delta_bias_b1, col="blue")
legend("topright", legend=c("b0", "b1"), col=c("green", "blue"), lty=1, cex=0.8)
dev.off()


b0_adj_var = b0_adj_sd^2
b1_adj_var = b1_adj_sd^2
b0_naive_var = b0_naive_sd^2
b1_naive_var = b1_naive_sd^2
delta_var_b0 = (b0_naive_var-b0_adj_var)/b0_naive_var
delta_var_b1 = (b1_naive_var-b1_adj_var)/b1_naive_var

min7 = min(delta_var_b0, delta_var_b1)
max7 = max(delta_var_b0, delta_var_b1)
dev.new()
pdf("delta_var_real_strc_size.pdf", 4, 4)
plot(graph_size_list, delta_var_b0, col="green", ylim=c(min7, max7), main = "(Var(Naive)-Var(Adj))/Var(Naive)", xlab="", ylab="")
lines(graph_size_list, delta_var_b0, col="green")
grid()
par(new=T)
plot(graph_size_list, delta_var_b1, col="blue", ylim=c(min7, max7), main = "", xlab="", ylab="")
lines(graph_size_list, delta_var_b1, col="blue")
legend("topright", legend=c("b0", "b1"), col=c("green", "blue"), lty=1, cex=0.8)
dev.off()

min8 = min(b1_adj_var, b1_naive_var)
max8 = max(b1_adj_var, b1_naive_var)
dev.new()
pdf("b1_real_strc_var_size.pdf", 4, 4)
plot(graph_size_list, b1_adj_var, col="green", ylim=c(min8, max8), main = "var(b1)", xlab="", ylab="")
lines(graph_size_list, b1_adj_var, col="green")
grid()
par(new=T)
plot(graph_size_list, b1_naive_var, col="red", ylim=c(min8, max8), main = "", xlab="", ylab="")
lines(graph_size_list, b1_naive_var, col="red")
legend("topright", legend=c("b1_adj", "b1_naive"), col=c("green", "red"), lty=1, cex=0.8)
dev.off()

dev.new()
pdf("Xhat1_bias_size.pdf", 4, 4)
plot(graph_size_list, mean_abs_X_bias_graphs, col="blue", main = "mean(abs(Xhat1-X1))", xlab="", ylab="")
lines(graph_size_list, mean_abs_X_bias_graphs, col="blue")
grid()
dev.off()


min3 = min(b1_true_sd, b1_adj_sd, b1_naive_sd)
max3 = max(b1_true_sd, b1_adj_sd, b1_naive_sd)
dev.new()
pdf("b1_real_strc_sd_size.pdf", 4, 4)
plot(graph_size_list, b1_true_sd, col="blue", ylim=c(min3, max3), main = "sd(b1)", xlab="", ylab="")
lines(graph_size_list, b1_true_sd, col="blue")
grid()
par(new=T)
plot(graph_size_list, b1_adj_sd, col="green", ylim=c(min3, max3), main = "", xlab="", ylab="")
lines(graph_size_list, b1_adj_sd, col="green")
par(new=T)
plot(graph_size_list, b1_naive_sd, col="red", ylim=c(min3, max3), main = "", xlab="", ylab="")
lines(graph_size_list, b1_naive_sd, col="red")
dev.off()

min4 = min(b0_true_sd, b0_adj_sd, b0_naive_sd)
max4 = max(b0_true_sd, b0_adj_sd, b0_naive_sd)
dev.new()
pdf("b0_real_strc_sd_size.pdf", 4, 4)
plot(graph_size_list, b0_true_sd, col="blue", ylim=c(min4, max4), main = "sd(b1)", xlab="", ylab="")
lines(graph_size_list, b0_true_sd, col="blue")
grid()
par(new=T)
plot(graph_size_list, b0_adj_sd, col="green", ylim=c(min4, max4), main = "", xlab="", ylab="")
lines(graph_size_list, b0_adj_sd, col="green")
par(new=T)
plot(graph_size_list, b0_naive_sd, col="red", ylim=c(min4, max4), main = "", xlab="", ylab="")
lines(graph_size_list, b0_naive_sd, col="red")
dev.off()


###########
###########
###########
###########
###########
###########
###########
###########
###########
#setwd("~/Desktop/redata/Feb6")
#dev.new()
#pdf("b0_estimate_rdpg_sim.pdf", 7, 5)
#par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b0_true, b0_naive,  b0_adj,  notch=TRUE, 
        main=bquote(paste("b0 Estimate","\n graph_size:", graph_size , "\n mc_runs:",length(b1_true))), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)
#dev.off()

#dev.new()
#pdf("b1_estimate_rdpg_sim.pdf", 7, 5)
boxplot(b1_true, b1_naive,  b1_adj,  notch=TRUE, 
        main=bquote(paste("b1 Estimate","\n graph_size:", graph_size , "\n mc_runs:",length(b1_true))), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)
#dev.off()


# df=data.frame(rbind(X, Xr_noise,  Xr_noise_simp, Xr_noise_hw ))
# df$type = c(rep(1, graph_size), rep(2, graph_size),  rep(3, graph_size), rep(4, graph_size))
# with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
# triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")


#dev.new()
#pdf("b0_MSE_rdpg_sim.pdf", 7, 5)
boxplot(b0_true_se, b0_naive_se,  b0_adj_se, notch=TRUE, 
        #   main=bquote(paste("b1 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)
#dev.off()

#dev.new()
#pdf("b1_MSE_rdpg_sim.pdf", 7, 5)
boxplot(b1_true_se, b1_naive_se,  b1_adj_se, notch=TRUE, 
        #   main=bquote(paste("b1 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)
#dev.off()



#CI 
t.test(b0_true)$conf.int
t.test(b0_naive)$conf.int
t.test(b0_adj)$conf.int

t.test(b1_true)$conf.int
t.test(b1_naive)$conf.int
t.test(b1_adj)$conf.int


#SE (sign test )
N = sum((b0_adj_se - b0_naive_se) !=0)
k = sum(b0_adj_se < b0_naive_se)
1 - pbinom(k, size=N, prob=0.5)

N = sum((b1_adj_se - b1_naive_se) !=0)
k = sum(b1_adj_se < b1_naive_se)
1 - pbinom(k, size=N, prob=0.5)








#generate an orthogonal matrix and rotate the X_hat by this matrix (Q)
#Q=(1/3)*matrix(c(2, 1, 2, -2, 2, 1, 1, 2, -2), nrow=3)
Q=diag(3)
Xr_noise = X_noise%*%Q

#rotate Xr_noise and project onto HW curve: 
Xr= X%*%Q
# 
# #Take Xr_noise onto HW, call this, X_naive 
# #Construct Q_est
# g=Xr_noise
# g_dist = dist(g, diag=T, upper=T)
# g_dist_m = as.matrix(g_dist, nrow=graph_size)
# which(g_dist_m == max(g_dist_m), arr.ind = TRUE)
# end1_index = which(g_dist_m == max(g_dist_m), arr.ind = TRUE)[1]
# end2_index = which(g_dist_m == max(g_dist_m), arr.ind = TRUE)[2]
# foo = g_dist_m[end1_index,]-g_dist_m[end2_index,]
# mid_index_z = which(foo==min(abs(foo)))
# if (is.integer(mid_index_z) && length(mid_index_z)==0) {mid_index_z = which(foo==-1*min(abs(foo)))[[1]]}
# 
# Xr_noise_pts = rbind(Xr_noise[end1_index,], Xr_noise[mid_index_z,], Xr_noise[end2_index,])
# 
# X_pts = rbind(X[1,], X[2,], X[3,])
# Xr_pts = rbind(Xr[1,], Xr[2,], Xr[3,])
# 
# Q_transpose_est_noise = t(solve(t(X_pts)%*%X_pts)%*%t(X_pts)%*%Xr_noise_pts)

# #Rotate X_hat back to simplex
# Xr_noise_simp = Xr_noise%*%Q_transpose_est_noise
# 
# #Project onto H-W
# min.RSS <- function(data, par) {
#   with(data, sum(((par^2-y[1])^2 + (2*par*(1-par)-y[2])^2 + ((1-par)^2-y[3])^2)^2))}
# 
# theta_noise_hw <- c()
# for (i in 1:nrow(X)){
#   dat=data.frame(y=Xr_noise_simp [i,])
#   result <- optim(par = 0.5, min.RSS, data = dat, method = "L-BFGS-B", lower = 0, upper = 1.0)
#   theta_noise_hw[i] <- result$par
# }
# 
# #CHEAT:
# flip = matrix(c(0, 0, 1, 0, 1, 0, 1, 0, 0), byrow=T, nrow=3)
# X_naive = theta_noise_hw
# if(sum((theta_noise_hw-theta_true)^2)>sum(((1-theta_noise_hw)-theta_true)^2)){
#   X_naive = 1-theta_noise_hw
#   theta_noise_hw = 1-theta_noise_hw
#   print("1")
#   Xr_noise_simp = (Xr_noise_simp) %*% flip
# }
# 

# #Construct estimate for covaraince (X_hat-X)
# sigma_u_theta_est_hw = c()
# for (i in 1:graph_size){
#   f = theta_noise_hw[i]
#   S11 = (1)*((1/126) + (f/15)+(f^2)/42-(f^3)/105-2*(f^4)/35)
#   S12 = (1)*((13/1260) + (4*f/105)-((f^2)/30)+(2*(f^3)/105)-((f^4)/70))
#   S13 = (1)*((1/180)+(f/105)-((f^2)/70)+(f^3)/105-(f^4)/210)
#   S22 = (1)*((1/45)+(4*f/105)-(2*f^2/35)+(4*f^3/105)-(2*f^4/105))
#   S23 = (1)*((5/252)+(f/35)-13*(f^2)/210+(4*f^3)/105-(f^4/70)) 
#   S33 = (1)*((2/63)+(f/7)-(73*(f^2)/210)+(5*(f^3)/21) - (2*(f^4)/35))
#   S_matrix = matrix(c(S11, S12, S13, S12, S22, S23, S13, S23, S33), nrow=3, byrow=T)
#   if(any(diag(S_matrix)<0)){print("ALERT")}
#   Sigma_X_est_hw = (1/graph_size)*Delta_inv%*%S_matrix%*%Delta_inv
#   sigma_u_theta_est_hw[i] = (1/4)*(1*Sigma_X_est_hw[1,1]-2*Sigma_X_est_hw[1,3]+Sigma_X_est_hw[3,3])
# }


