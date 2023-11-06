library(igraph)
library(netrankr)
library(network)
library(graphkernels)
library(arrangements)

####Generating an ERMG####
sample_ermg = function(n1,n2,p){
  # Sampling within groups
  gw1 = sample_gnp(n1, p[1], directed = FALSE, loops = FALSE)
  #plot(gw1)
  n_e11=gsize(gw1)
  E(gw1)$group <- rep(1,n_e11)
  #get.edge.attribute(gw1)
  
  gw2 = sample_gnp(n2, p[3], directed = FALSE, loops = FALSE)
  #plot(gw2)
  n_e22=gsize(gw2)
  E(gw2)$group <- rep(2,n_e22)
  #get.edge.attribute(gw2)
  
  #Between two groups
  gb= sample_bipartite(n1, n2, p = p[2])
  #plot(gb)
  n_e12=gsize(gb)
  E(gb)$group <- rep(12,n_e12)
  #get.edge.attribute(gb)
  
  gw = gw1 + gw2
  #plot(gw)
  gg <- gb %u% gw
  #plot(gg)
  #get.edge.attribute(gg)
  
  E(gg)$group_1[is.na(E(gg)$group_1)] <- 0
  E(gg)$group_2[is.na(E(gg)$group_2)] <- 0
  E(gg)$group <- E(gg)$group_1 + E(gg)$group_2
  #get.edge.attribute(gg)
  
  G = gg
  n=n1+n2
  V(G)$name <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[1:n]
  G
}

########


sample_ermg_with_ecount = function(n1,n2,p){
  # Sampling within groups
  gw1 = sample_gnp(n1, p[1], directed = FALSE, loops = FALSE)
  #plot(gw1)
  n_e11=gsize(gw1)
  E(gw1)$group <- rep(1,n_e11)
  #get.edge.attribute(gw1)
  
  gw2 = sample_gnp(n2, p[3], directed = FALSE, loops = FALSE)
  #plot(gw2)
  n_e22=gsize(gw2)
  E(gw2)$group <- rep(2,n_e22)
  #get.edge.attribute(gw2)
  
  #Between two groups
  gb= sample_bipartite(n1, n2, p = p[2])
  #plot(gb)
  n_e12=gsize(gb)
  E(gb)$group <- rep(12,n_e12)
  #get.edge.attribute(gb)
  
  gw = gw1 + gw2
  #plot(gw)
  gg <- gb %u% gw
  #plot(gg)
  #get.edge.attribute(gg)
  
  E(gg)$group_1[is.na(E(gg)$group_1)] <- 0
  E(gg)$group_2[is.na(E(gg)$group_2)] <- 0
  E(gg)$group <- E(gg)$group_1 + E(gg)$group_2
  #get.edge.attribute(gg)
  
  G = gg
  n=n1+n2
  V(G)$name <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[1:n]
  list(G = G, n_e11 = n_e11, n_e12 = n_e12, n_e22 = n_e22 )
}

##########Sampling ERMG with cliques#################

sample_ermg_planted_clique = function(n1,n2,p,K){
  if(K==0 || K==1){ G = sample_ermg(n1,n2,p)}
  else
    if(K>1){
  G = sample_ermg(n1,n2,p)
  smpl = sample(vertex_attr(G)$name,K)
  new_edges = as.vector(combinations(smpl, k = 2, layout = "column"))
  G <- G + edges(new_edges)
  l = get.edge.ids(G, new_edges)
  t = vertex_attr(G, index = new_edges)$type
  q = seq(1,2*length(l),2)
  for (k in 1:length(l)) {
    
    if(t[q[k]] == FALSE && t[q[k]+1] == FALSE)E(G)$group[l[k]] = 1
    else
      if (t[q[k]] == TRUE && t[q[k]+1] == FALSE)E(G)$group[l[k]] = 12
      else
        if(t[q[k]] == FALSE && t[q[k]+1] == TRUE)E(G)$group[l[k]] = 12
        else
          if(t[q[k]] == TRUE && t[q[k]+1] == TRUE)E(G)$group[l[k]] = 2
  }
  G = simplify(G, remove.multiple = TRUE, edge.attr.comb="first")
    }
} 

##############

sample_ermg_planted_cliques = function(n1,n2,p,K,M){
  if(K==0 || K==1){ G = sample_ermg(n1,n2,p)}
  else
    if(K>1){
      G = sample_ermg(n1,n2,p)
      gg = G
      for(i in 1:M){
        smpl = sample(vertex_attr(G)$name,K)
        new_edges = as.vector(combinations(smpl, k = 2, layout = "column"))
        gg <- gg + edges(new_edges)
        l = get.edge.ids(gg, new_edges)
        t = vertex_attr(gg, index = new_edges)$type
        q = seq(1,2*length(l),2)
        for (k in 1:length(l)) {
          if(t[q[k]] == FALSE && t[q[k]+1] == FALSE)E(gg)$group[l[k]] = 1
          else
            if (t[q[k]] == TRUE && t[q[k]+1] == FALSE)E(gg)$group[l[k]] = 12
            else
              if(t[q[k]] == FALSE && t[q[k]+1] == TRUE)E(gg)$group[l[k]] = 12
              else
                if(t[q[k]] == TRUE && t[q[k]+1] == TRUE)E(gg)$group[l[k]] = 2
        }
        gg = simplify(gg, remove.multiple = TRUE, edge.attr.comb="first")
      }
    }
  list(Original = G, with_cliques = gg)
}  

##############Edge group assignment#########

edge_group_assign = function(G,v_group){
  r = seq(1,2*gsize(G),2)
  n = 0
  for(i in 1:gsize(G)){
    v = as.vector(V(G)[.inc(i)])
    n[r[i]:(r[i]+1)] = v
  }
  V(G)$group = v_group
  E(G)$group = 0
  l= get.edge.ids(G, n)
  t = vertex_attr(G, index = n)$group
  q = seq(1,2*length(l),2)
  for (k in 1:length(l)) {
    if(t[q[k]] == 1 && t[q[k]+1] == 1)E(G)$group[l[k]] = 1
    else
      if (t[q[k]] == 1 && t[q[k]+1] == 2)E(G)$group[l[k]] = 12
      else
        if(t[q[k]] == 2 && t[q[k]+1] == 1)E(G)$group[l[k]] = 12
        else
          if(t[q[k]] == 2 && t[q[k]+1] == 2)E(G)$group[l[k]] = 2
  }
  G
}

################Edge counter#########
edge_counter_ermg = function(group)
{
  n_e11 = 0
  n_e12 = 0
  n_e22 = 0
  
  for (i in 1:length(group)) {
    if(group[i] == 1) n_e11 = n_e11 + 1
    else
      if(group[i] == 12) n_e12 = n_e12 + 1
      else
        if(group[i] == 2) n_e22 = n_e22 + 1
  }
  list(n_e11 = n_e11, n_e12 = n_e12, n_e22 = n_e22 )
}


###############t_fun####################

p=c()
  t_fun = function(G){
    X = as_adj(G, edges = TRUE)
    p_ij_matrix = matrix(0, nrow(X), ncol(X))
    
    for (i in 1:ncol(X)) {
      for (j in 1:nrow(X)) {
        x = X[i,j]
        if(x ==0) p_ij_matrix[i,j] = 0
        else
          if(E(G)[x]$group == 1) p_ij_matrix[i,j] = p[1]
          else
            if(E(G)[x]$group == 2) p_ij_matrix[i,j] = p[3]
            else
              if(E(G)[x]$group == 12) p_ij_matrix[i,j] = p[2]
      }
    }
    p_ij_matrix
  }
  
  
  ##Kernel Computation####
  compute.transition.list = function(X){
    P=list()
    for (w in 1:length(X)){
      x = X
      x[w] = abs(1-X[w])
      G = graph_from_adjacency_matrix(x)
      P[[w]] = G
    }
    P[[length(X)+1]] = graph_from_adjacency_matrix(X)
    P
  }
  
  
  compute.sampled.list = function(X, sample.index){
    P=list()
    l = length(sample.index)
    for (w in 1:l){
      x = X
      x[sample.index[w]] = abs(1 - X[sample.index[w]])
      G = graph_from_adjacency_matrix(x)
      P[[w]] = G
    }
    P[[l+1]] = graph_from_adjacency_matrix(X)
    P
  }
  
  
  compute.normalized = function(K){
    V = diag(K)
    D.mat = diag(1.0/sqrt(V))
    D.mat%*%K%*%D.mat
  }
  
  
  ## Weisfeiler_Lehman Graph Kernel
  compute.wl.kernel=function(X, level=3, diag=1, normalize=TRUE){
    n = length(X)
    P = compute.transition.list(X)
    kernel.matrix = CalculateWLKernel(P, level)
    rm(P)
    K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
    K.vec = kernel.matrix[1:n,n+1]
    K = K + outer(K.vec, K.vec, function(x,y)x+y)
    if(normalize)K = compute.normalized(K)
    if(diag==0)diag(K)<-0
    K
  }
  
  # GKSS methods
  
  
  # Using vvRKHS
  generate.one.GKSS.condition=function(G, X, kernel=compute.wl.kernel,diagonal=1)
  {
    S=matrix(X,byrow=TRUE)
    n.sim=length(S)
    S.t =  t_fun(G)
    S.t.vec=abs(S - matrix(S.t, byrow=TRUE))
    S.mat = S.t.vec %*% t(S.t.vec)
    Ky.mat = (S*2-1)%*%t(S*2-1) 
    
    
    K = kernel(X, diag=diagonal)
    J.kernel = S.mat * Ky.mat * K
    
    W=rep(1/n.sim,n.sim)
    J.kernel.out=J.kernel
    if(diagonal==0)diag(J.kernel)=rep(0,n.sim)
    v.kernel = var(S)
    stats.value=n.sim*t(W)%*%J.kernel%*%W * sqrt(v.kernel)
    #Return:
    #stats.value: n times KSD
    #J.kernel: J kernel matrix for wild bootstrapt 
    list(stats.value=stats.value,J.kernel=J.kernel.out, K=K, S=S.mat)
  }

##################chi_square_BC####################
  
  ##v_group should have values 1 or 2
  ##m_rl = c(m11,m12,m22)
  ##n_l = c(n1,n2)
  chi_square_BC = function(G,q,n_l,m_rl,v_group)
  {
    n1= n_l[1]
    n2= n_l[2]
    n = n1 + n2
    V(G)$group = v_group
    
    N = matrix(0,nrow = n, ncol = q)
    for (i in 1:n) {
      nei = neighbors(G, i)
      if(length(nei)==0) N[i,]=0
      else
        if(length(nei) > 0){
          for (j in 1:length(nei)) {
            if(nei[j]$group == 1) N[i,1] = N[i,1]+1
            else
              if(nei[j]$group == 2) N[i,2] = N[i,2]+1
          }
        }
    }
    
    m11 = m_rl[1]
    m12 = m_rl[2]
    m22 = m_rl[3]
    v1 = c(2*(m11/n1), m12/n1)
    v2 = c(m12/n2, 2*(m22/n2))
    
    E = matrix(0,nrow = n, ncol = q)
    for (i in 1:n) {
      if(V(G)[i]$group == 1) E[i,] = v1
      if(V(G)[i]$group == 2) E[i,] = v2
    }
    
    my_fun = function(x,y) {((x-y)^2)/y}
    ts = my_fun(N,E)
    test_stat = sum(ts)
    test_stat
  }