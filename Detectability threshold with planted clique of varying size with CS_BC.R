library(igraph)
library(netrankr)
library(network)
library(graphkernels)
library(arrangements)

alpha=0.05
n = 20
n1=10
n2=10


cout = 50
cin_1 =  1 + cout - sqrt(1 + 4 *cout)
cin_2 = 1 + cout + sqrt(1 + 4 *cout)
n_c=100

cin = 45

p_0 = c(cin/n_c,cout/n_c,cin/n_c)

N = 200
m = 100

################################
######Null Model Parameters#####
################################

#K = 0

##################################
###Alternative Model Parameters###
##################################

K = seq(0,20,1)


  ##############################
  ########Methodology###########
  ##############################
  
  # re-direct to the desired directory
  setwd("C:/Users/fatima/OneDrive - Nexus365/DPhil/Research/Code books/Networks") #for windows working directory
  source("Methodology_for_simulations_code.R")

#######################################################
##############Simulations from null model##############
#######################################################
{
  statistic = c()
  p = p_0
  
  for (i in 1:N){
    G = sample_ermg(n1,n2,p)
    X = G[,]
    #t_fun(G)
    statistic[i] = generate.one.GKSS.condition(G, X, kernel=compute.wl.kernel,diagonal=1)$stats.value
  }
  
  critical_value_lower = quantile(statistic,alpha/2)
  critical_value_upper = quantile(statistic,(1-alpha/2))
  mean(statistic)
}

###############Chi_square_BC########################
n_l = c(n1,n1)
p = p_0

test_stat_CS_BC = c()
for (i in 1:N){
  g = sample_ermg(n_l[1],n_l[2],p)
  E = edge_counter_ermg(E(g)$group)
  m_rl = c(E$n_e11,E$n_e12,E$n_e22)
  
  v_group = replace(V(g)$type, V(g)$type == "TRUE", 2 )
  v_group = replace(v_group, v_group == 0, 1 )
  
  #test_stat_CS_BC[i] = 
  test_stat_CS_BC[i] = chi_square_BC(g,2,n_l,m_rl,v_group)
}

#test_stat_CS_BC
critical_value_CS_BC_lower = quantile(test_stat_CS_BC,alpha/2)
critical_value_CS_BC_upper = quantile(test_stat_CS_BC,(1-alpha/2))
mean(test_stat_CS_BC)

####################################################
##Computing test statistic under alternative model##
####################################################

p = p_0
l= length(K)
power_gKSS = matrix(0,l)
power_GLRT = matrix(0,l)
power_CSBC = matrix(0,l)

for (i in 1:l) {
  
  p_1 = p
  
  #######################################################
  test_stat.list_gKSS = matrix(0,m)
  decision_gKSS = matrix(0,m)
  test_stat.list_GLRT = matrix(0,m)
  decision_GLRT = matrix(0,m)
  test_stat.list_CSBC = matrix(0,m)
  decision_CSBC = matrix(0,m)
  
  for (j in 1:m) {
    
    test_G = sample_ermg_planted_clique(n1,n2,p_1,K[i])
    Y = test_G[,]
    
    ########gKSS###########
    
    test_stat.list_gKSS[j] = generate.one.GKSS.condition(test_G, Y, kernel=compute.wl.kernel,diagonal=1)$stats.value
    
    decision_gKSS[j] = isTRUE(critical_value_lower >= test_stat.list_gKSS[j] ||test_stat.list_gKSS[j] >= critical_value_upper)
    #False means do not reject null
    
    #########GLRT##########
    
    ############MLE Paramter estimates############
    
    E = edge_counter_ermg(E(test_G)$group)
    a_e11 = E$n_e11
    a_e12 = E$n_e12
    a_e22 = E$n_e22
    
    #p1_h = a_e11/(n1*(n1-1)/2)
    #p12_h = a_e12/ (n1*n2)
    #p2_h = a_e22/ (n2*(n2-1)/2)
    p_h = (a_e11+a_e22)/((n1*(n1-1)/2)+(n2*(n2-1)/2))
    
    #############Test Statistic###################    
    
    L_null = ((p[1])^(a_e11))*((1-p[1])^(n1*(n1-1)/2 - a_e11)) * ((p[2])^(a_e12))*((1-p[2])^(n1*n2 - a_e12)) * ((p[3])^(a_e22))*((1-p[3])^(n2*(n2-1)/2 - a_e22))
    
    L_alt = ((p_h)^(a_e11))*((1-p_h)^(n1*(n1-1)/2 - a_e11)) * ((p[2])^(a_e12))*((1-p[2])^(n1*n2 - a_e12)) * ((p_h)^(a_e22))*((1-p_h)^(n2*(n2-1)/2 - a_e22))
    
    GLR = L_null/L_alt
    
    test_stat.list_GLRT[j] = -2*log(GLR)
    
    decision_GLRT[j] = isTRUE(test_stat.list_GLRT[j] > qchisq(1 - alpha, df=1) )
    #False means do not reject null  
    
    ##########CSBC####################
    m_rl = c(E$n_e11,E$n_e12,E$n_e22)
    v_group = replace(V(test_G)$type, V(test_G)$type == "TRUE", 2 )
    v_group = replace(v_group, v_group == 0, 1 )
    
    test_stat.list_CSBC[j] =chi_square_BC(test_G,2,n_l,m_rl,v_group)
    
    decision_CSBC[j] = isTRUE(critical_value_CS_BC_lower >= test_stat.list_CSBC[j] ||test_stat.list_CSBC[j] >= critical_value_CS_BC_upper)
    
  }
  power_gKSS[i] = mean(decision_gKSS)
  power_GLRT[i] = mean(decision_GLRT)
  power_CSBC[i] = mean(decision_CSBC)
}
power_gKSS_3 = power_gKSS
power_GLRT_3 = power_GLRT
power_CSBC_3 = power_CSBC

#####################Plot############################
plot(K,power_gKSS_3, xlab = "K", ylab = "Rejection Rate",ylim = c(0,1),pch=19)
lines(K,power_gKSS_3, lty =3, col = 1)
#abline(v= cout/n_c, col = "blue" )
#abline(v= cin_1/n_c, col = "red")
#abline(v= cin_2/n_c, col = "red")
points(K, power_GLRT_3, lty = 3, col = 2, pch=19)
lines(K, power_GLRT_3, lty = 3, col = 2)
points(K, power_CSBC_3, lty = 3, col = "blue", pch=19)
lines(K, power_CSBC_3, lty = 3, col = "blue")

abline(h = alpha,col = "gold")

legend("topleft", inset=0.01,title = "Test",legend= c("gKSS","GLRT", "Chi_Square_BC"), lty = 3, col = c(1,2,"blue") ,bty = "o",cex=0.9, pt.cex = 1)

#########################################################################

# power_gKSS_2, power_GLRT_2 takes c_out = 50, c_in = 30
# power_gKSS_3, power_GLRT_3 takes c_out = 50, c_in = 45
# power_gKSS_4, power_GLRT_4 takes c_out = 50, c_in = 60
# power_gKSS_5, power_GLRT_5 takes c_out = 50, c_in = 80

