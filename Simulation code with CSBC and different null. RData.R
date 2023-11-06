library(igraph)
library(netrankr)
library(network)
library(graphkernels)

alpha=0.05
n = 20
n1=10 
n2=10

N = 200
m = 100

################################
######Null Model Parameters#####
################################

#p_0 = c(p1,p12,p2)
p_0 = c(0.5,0.5,0.5)

##################################
###Alternative Model Parameters###
##################################

p_alt = p_0
p12_per = seq(0, 1, 0.05)

##############################
########Methodology###########
##############################

# re-direct to the desired directory
setwd("C:/Users/fatima/OneDrive - Nexus365/DPhil/Research/Code books/Networks") #for windows working directory
source("Methodology_for_simulations_code.R")

#######################################################
##############Simulations from null model##############
#######################################################

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
l= length(p12_per)
power_gKSS = matrix(0,l)
power_GLRT = matrix(0,l)
power_CSBC = matrix(0,l)

for (i in 1:l) {
  
  p_1 = c(p_alt[1],p12_per[i],p_alt[3])
  
  #######################################################
  test_stat.list_gKSS = matrix(0,m)
  decision_gKSS = matrix(0,m)
  test_stat.list_GLRT = matrix(0,m)
  decision_GLRT = matrix(0,m)
  test_stat.list_CSBC = matrix(0,m)
  decision_CSBC = matrix(0,m)
  
  for (j in 1:m) {
    
    test_G = sample_ermg(n1,n2,p_1)
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
      p12_h = a_e12/ (n1*n2)
      #p2_h = a_e22/ (n2*(n2-1)/2)
    
    #############Test Statistic###################    
    
    L_null = ((p[1])^(a_e11))*((1-p[1])^(n1*(n1-1)/2 - a_e11)) * ((p[2])^(a_e12))*((1-p[2])^(n1*n2 - a_e12)) * ((p[3])^(a_e22))*((1-p[3])^(n2*(n2-1)/2 - a_e22))
    
    L_alt = ((p[1])^(a_e11))*((1-p[1])^(n1*(n1-1)/2 - a_e11)) * ((p12_h)^(a_e12))*((1-p12_h)^(n1*n2 - a_e12)) * ((p[3])^(a_e22))*((1-p[3])^(n2*(n2-1)/2 - a_e22))
    
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
  #print(test_stat.list_CSBC)
  #print(decision_CSBC)
  power_gKSS[i] = mean(decision_gKSS)
  power_GLRT[i] = mean(decision_GLRT)
  power_CSBC[i] = mean(decision_CSBC)
}
power_gKSS_3 = power_gKSS
power_GLRT_3 = power_GLRT
power_CSBC_3 = power_CSBC
###################################################################
par(mar = c(3, 3, 2, 6),xpd=TRUE)
plot(p12_per,power_gKSS_1, xlab = "Between Groups Edge Probability", ylab = "Rejection Rate",type = "n", ylim = c(0,1))
lines(p12_per,power_gKSS_1, lty = 1)
lines(p12_per,power_GLRT_1, lty = 2)
#lines(p12_per,power_CSBC_1, lty = 3)

lines(p12_per,power_gKSS_2, lty = 1, col = 2)
lines(p12_per,power_GLRT_2, lty = 2, col = 2)
#lines(p12_per,power_CSBC_2, lty = 3, col = 2)

lines(p12_per,power_gKSS_3, lty = 1, col = 3)
lines(p12_per,power_GLRT_3, lty = 2, col = 3)
#lines(p12_per,power_CSBC_3, lty = 3, col = 3)

lines(p12_per,power_gKSS_4, lty = 1, col = 4)
lines(p12_per,power_GLRT_4, lty = 2, col = 4)
#lines(p12_per,power_CSBC_4, lty = 3, col = 4)

lines(p12_per,power_gKSS_5, lty = 1, col = 6)
lines(p12_per,power_GLRT_5, lty = 2, col = 6)
#lines(p12_per,power_CSBC_5, lty = 3, col = 6)

legend("right", inset = - 0.21 ,title = "Null Model",legend= c("(0.08, 0.1 , 0.09)","(0.1,0.3,0.2)","(0.5,0.5,0.5)","(0.4,0.7,0.6)","(0.9,0.9,0.6)") ,bty = "n",cex=0.7, pt.cex = 1, fill = c(seq(1,4,1),6))
legend("topright", inset = c(-0.16,  0.2) ,title = "Test",legend= c("gKSS", "GLRT") ,bty = "n",cex=0.7, pt.cex = 2, lty = c(1,2))

par(xpd=FALSE)
abline(h = alpha,col = "gold")
###########################################################################################

#Scenerio 1 : (0.08, 0.1 , 0.09)
#Scenerio 2 : (0.1,0.3,0.2)
#Scenerio 3 : (0.5,0.5,0.5)
#Scenerio 4 : (0.4,0.7,0.6)
#Scenerio 5 : (0.9,0.9,0.6)