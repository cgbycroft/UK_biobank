# This hmm function(s) is based on the one created by Clare in /well/donnelly/ukbiobank_project_8874/clare/sex/scripts/ClaresHMM2.R
# It is updated to infer copy number 1, 2, or 3.

# trans is required. E.g
# stayProb = 0.999
#trans=diag(stayProb,dim(stateProbabilities)[3])  # make it hard to change state
#trans[1,2:3] = (1-stayProb)*c(0.5,0.5)
#trans[2,c(1,3)] = (1-stayProb)*c(0.5,0.5)
#trans[3,c(1,2)] = (1-stayProb)*c(0.5,0.5)
#trans=rbind(c(0.999,0.001),c(0.001,0.999))


myForwardsBackwards.OLD <- function(log2Ratios,trans,trans_0=NULL,scale_trans=NULL,HMMPosteriors){
    print('start hmm')
    print(Sys.time())
# trans_0 = initial prob.
# trans = transition matrix (nStates x nStates)
# log2Ratios = array:  nSNPs x nSamples
# HMMPosteriors = array: nSNPs x nStates x 2
# NOTE: The SNPs in Log2Ratios and HMMPosteriors must be in the same order!    
# scale_a = an optional vector of length nObservations - 1 to scale the baseline transition matrix (a) by its value at each observation towards a 50/50 transition. 
  # Values are in [0,1], where 1 is keep the baseline transition matrix

    
    nObs = dim(log2Ratios)[1]
    nInds = dim(log2Ratios)[2]
    nStates = dim(trans)[2]
    
    print(nStates)
    
    # set defaults
    if(is.null(trans_0)) trans_0 = matrix(1/nStates,nStates,nStates)
    
    # first compute emissions for each state at each position (for use later). Note that posterior list should be in the same order as the observations
    print('computing emissions...')

    # emissions is a list: nSNPs x nSamples x 2. Compute over SNPs.
    emissions = sapply(1:nObs,FUN=function(p){
           if(p%%1000==0) {
                print(p)
                print(Sys.time())
            }
        getEmissions(x=log2Ratios[p,],HMMPosteriors[p,,])
    },simplify=F)

    
    states = 1:nStates
      
    #maxScale = 0.5/trans[1,] #maximum scaling
    #minScale = c(1,1) #minumum scaling
    maxScale = (1/nStates)/trans[1,] #maximum scaling
    minScale = rep(1,nStates) #minumum scaling

    # How to get from snp distance to scaled transition probabilities. The matrix trans is fixed across all SNPs. In theory this could be learned by Baum-Welch.
    scaleFunction = function(s) {
     s2 =  maxScale[2]-(maxScale[2] - minScale[2])*s
     s1 = (1 - trans[1,2]*s2)/trans[1,1]
     if(sum(trans[1,]*c(s1,s2))!=1) print("error")
     return(c(s1,s2))
    }

    
    # forward part of the algorithm
    print('computing forwards...')

    f_norm = matrix(NA,nrow=nInds,ncol=nObs)
    fwd = array(dim = c(nObs,nInds,length(states)))
    f_curr_empty = matrix(NA,nrow=nInds,ncol=nStates) 
    
    for(i in 1:nObs){
        # print(i)
        #x_i = obs[i,]
        #n = N[i,]
        emis_i = emissions[[i]]  # emis_i is a matrix nInds x nStates
        f_curr = f_curr_empty # matrix with rows nInds and columns number of states
        
        if(i > 1){
            if(!is.null(scale_trans)){              
              scale_t = scaleFunction(scale_trans[i-1]) # first value in scale relates to 2nd iteration 
              #print(paste('scaling transition matrix by ',scale_t))
              trans2 = rbind(scale_t*trans[1,],rev(scale_t)*trans[2,])
            } else {
              trans2 = trans
            }
            if(sum(rowSums(trans2)==1)!=nStates) print('error with trans2')
        }
          
        
        for (st in states){
            if (i == 1) {
              prev_f_sum = rep(trans_0[st],nInds) # base case for the forward part. vector length nInds
            } else {              
              prev_f_sum = f_prev%*%trans2[,st] # vector length nInds
            }
            
            f_curr[,st] = emis_i[,st] * prev_f_sum  # get likelihoods for this observation in state st
        }        
        
        if(sum(f_curr < 0)>0) print('error: negative likelihoods')
                
        f_curr_sum = rowSums(f_curr)  # vector of length nInds
        
        if(sum(1/f_curr_sum > 10000000)>0){
          print(i)
          print(prev_f_sum)
          print(f_curr)
          print('')
        }
        
        f_norm[,i] = 1/f_curr_sum #store normalising factors for later. matrix of dim nInds x nObs

        f_curr = f_curr/f_curr_sum  # normalising to make sum == 1  (matrix of nInds x nStates)
        
        fwd[i,,] = f_curr  # append current ai to previous one
        
        f_prev = f_curr
    }

    
   # backwards part of the algorithm
    print('computing backwards...')  # for some reason this is much slower...

    bkw = array(dim = c(nObs,nInds,nStates))
    scale_trans_rev = rev(scale_trans)
    b_curr_empty = matrix(NA,nrow=nInds,ncol=nStates)# matrix of nInds x number of states

    emissionsRev = emissions[nObs:1] # reverse the emissions so going backwards
    f_norm_rev = f_norm[,nObs:1]
    
    for( i in 1:nObs){
#        for( i in 1:50){
        if(i%%1000==0) {
            print(i)
            print(Sys.time())
        }
        
        emis_i_plus = emissionsRev[[i]]  # emission from end of observations
 
        b_curr = b_curr_empty # matrix of nInds x number of states
        
        if (i > 1){ 
            if(!is.null(scale_trans)) {
              scale_t = scaleFunction(scale_trans_rev[i-1])
              #print(paste('scaling by',scale_t))
              trans2 = rbind(scale_t*trans[1,],rev(scale_t)*trans[2,])              
            } else {
              trans2 = trans
            }
            if(sum(rowSums(trans2)==1)!=nStates) print('error with trans2')
        }
                       
        for (st in states){
            if (i == 1) {
              b_curr[,st] = rep(1,nInds)
            } else {
              e=(emis_i_plus * b_prev)%*%trans2[st, ] # sum over states, after multiplying emissions by previous values              
              b_curr[,st] = e
            }
        }

#        b_curr = b_curr*f_norm[,nObs:1][,i] # %%% normalising by same scale as forwards part
         b_curr = b_curr*f_norm_rev[,i] # normalising by same scale as forwards part
        
        bkw[(nObs - i + 1),,] = b_curr # append bi on to the front of bkw
        
        b_prev = b_curr
    }

    print('computing snp-wise posteriors...')
        
    p_fwd_log = -rowSums(log10(f_norm))

    # get posterior probability of being in some state at some snp, given all the data.
    posterior = sapply(1:nInds,function(i){
        fwd[,i,]*bkw[,i,]/rowSums(fwd[,i,]*bkw[,i,])

    },simplify=F)

    print('return hmm')
    print(Sys.time())  

    return(posterior)
}



# This version allows updates of transition probabilities, and computes these more efficiently.
#HMMPosteriors = hmmResultsRound$parameters[[1]]

getTransitionsAndEmissions <- function(log2Ratios,trans,scale_trans=NULL,HMMPosteriors){
# this is required for viterbi and forwards backwards
    
    nObs = dim(log2Ratios)[1]
    nInds = dim(log2Ratios)[2]
    nStates = dim(trans)[2]
    
    print(nStates)
        
    # first compute emissions for each state at each position (for use later). Note that posterior list should be in the same order as the observations
    print('computing emissions...')

    # emissions is an array: nSamples x nStates x nSNPs. Computed over SNPs.
    emissions = sapply(1:nObs,FUN=function(p){
           if(p%%1000==0) {
                print(p)
                print(Sys.time())
            }
        getEmissions(x=log2Ratios[p,],HMMPosteriors[p,,])
    },simplify="array")

    
    # How to get from snp distance to scaled transition probabilities. The matrix trans is fixed across all SNPs. In theory this could be learned by Baum-Welch.

    scaleFunction = function(s){
        stateTransitions = sapply(1:nStates,function(l){
            ps = trans[l,]
            sapply(1:nStates,function(j){
                if(l==j) prob = 1 - sum(ps[-j]*( 1 - s ))
                if(l!=j) prob = ps[j]*( 1 - s )
                return(prob)
            },simplify="array")                   
        },simplify="array")
        return(t(stateTransitions))
    }

    # apply this across all snps
    if(is.null(scale_trans)) scale_trans = rep( 0.5,nObs-1 ) # assume 0.5 exponential distance between sites.
    transitions = sapply(scale_trans,scaleFunction,simplify="array")
    # transitions is array: nStates x nStates x nObs-1 (one matrix for each gap between snps)

    return(list("transitions"=transitions,"emissions"=emissions))

}


myForwardsBackwards <- function(log2Ratios,trans,trans_0=NULL,scale_trans=NULL,HMMPosteriors,computeTransitions=TRUE,computeViterbi=FALSE){
    
    print('start hmm - forwards-backwards')
    print(Sys.time())
# trans_0 = initial prob.
# trans = transition matrix (nStates x nStates)
# log2Ratios = array:  nSNPs x nSamples
# HMMPosteriors = array: nSNPs x nStates x 2
# NOTE: The SNPs in Log2Ratios and HMMPosteriors must be in the same order!    
# scale_a = an optional vector of length nObservations - 1 to scale the baseline transition matrix (a) by its value at each observation towards a 50/50 transition. 
  # Values are in [0,1], where 1 is keep the baseline transition matrix

    
    nObs = dim(log2Ratios)[1]
    nInds = dim(log2Ratios)[2]
    nStates = dim(trans)[2]
    states = 1:nStates
    
    print(nStates)

    # get transitions and emissions
    TE = getTransitionsAndEmissions(log2Ratios,trans,scale_trans,HMMPosteriors)
    
    transitions = TE$transitions
    emissions = TE$emissions        

    # set default for initial trans
    if(is.null(trans_0)) trans_0 = matrix(1/nStates,nStates,nStates)
    

    # forward part of the algorithm
    print('computing forwards...')

    f_norm = matrix(NA,nrow=nInds,ncol=nObs)
    fwd = array(dim = c(nObs,nInds,length(states)))
    f_curr_empty = matrix(NA,nrow=nInds,ncol=nStates) 
    
    for(i in 1:nObs){
        # print(i)
        #x_i = obs[i,]
        #n = N[i,]
        emis_i = emissions[,,i]  # emis_i is a matrix nInds x nStates
        f_curr = f_curr_empty # matrix with rows nInds and columns number of states
        
        if(i > 1){
            # get transition probabilities
            trans2 = transitions[,,(i-1)]
            if(sum(rowSums(trans2)==1)!=nStates) print( paste0('error with trans2 ',i) )
        }
          
        
        for (st in states){
            if (i == 1) {
              prev_f_sum = rep(trans_0[st],nInds) # base case for the forward part. vector length nInds
            } else {              
              prev_f_sum = f_prev%*%trans2[,st] # vector length nInds
            }
            
            f_curr[,st] = emis_i[,st] * prev_f_sum  # get likelihoods for this observation in state st
        }        
        
        if(sum(f_curr < 0)>0) print( paste0('error: negative likelihoods at observation ',i) )
                
        f_curr_sum = rowSums(f_curr)  # vector of length nInds
        
        if(sum(1/f_curr_sum > 10000000)>0){
          print(i)
          print(prev_f_sum)
          print(f_curr)
          print('')
        }
        
        f_norm[,i] = 1/f_curr_sum #store normalising factors for later. matrix of dim nInds x nObs

        f_curr = f_curr/f_curr_sum  # normalising to make sum == 1  (matrix of nInds x nStates)
        
        fwd[i,,] = f_curr  # append current ai to previous one
        
        f_prev = f_curr
    }

    
   # backwards part of the algorithm
    print('computing backwards...')  # for some reason this is much slower...

    bkw = array(dim = c(nObs,nInds,nStates))
    transitions_rev = transitions[,,(nObs-1):1]
    b_curr_empty = matrix(NA,nrow=nInds,ncol=nStates)# matrix of nInds x number of states

    emissions_rev = emissions[,,nObs:1] # reverse the emissions so going backwards
    f_norm_rev = f_norm[,nObs:1]
    
    for( i in 1:nObs){
#        for( i in 1:50){
        if(i%%1000==0) {
            print(i)
            print(Sys.time())
        }
        
        emis_i_plus = emissions_rev[,,i]  # emission from end of observations
 
        b_curr = b_curr_empty # matrix of nInds x number of states
        
        if (i > 1){
            # get transition probabilities (but reverse). NOTE: assume forwards/backwards symmetry in transitions
            trans2 = transitions_rev[,,(i-1)]
            if(sum(rowSums(trans2)==1)!=nStates) print( paste0('error with trans2 ',i)  )
        }
                       
        for (st in states){
            if (i == 1) {
              b_curr[,st] = rep(1,nInds)
            } else {
              e=(emis_i_plus * b_prev)%*%trans2[st, ] # sum over states, after multiplying emissions by previous values              
              b_curr[,st] = e
            }
        }

         b_curr = b_curr*f_norm_rev[,i] # normalising by same scale as forwards part
        
        bkw[(nObs - i + 1),,] = b_curr # append bi on to the front of bkw
        
        b_prev = b_curr
    }

    
    print('computing snp-wise posteriors...')
        
    p_fwd_log = -rowSums(log10(f_norm))

    # get posterior probability of being in some state at some snp, given all the data.
    posterior = sapply(1:nInds,function(k){
        fwd[,k,]*bkw[,k,]/rowSums(fwd[,k,]*bkw[,k,])

    },simplify="array")

    # get posterior probability of transition from one state to another at the next snp. Formula from stanford tutorial. Only required if need to update transition probabilities.
    transitionsPosterior = NULL
    if( computeTransitions ){

        print("Computing posterior transition probabilities for each individual...")
        
        distanceScale = 1/(1 - scale_trans)
        distanceScale = length(distanceScale) * distanceScale/sum(distanceScale)
            
        transitionsPosterior = sapply(1:nInds,function(k){
            mat = sapply(states,function(i){
                
                sapply(states,function(j){
                    
                    #postt = posterior[1:(nObs - 1),i,k]  # posteriors for ind k at state i at t
                    aij = transitions[i,j,] # transitions i to j (scaled by distance)
                    ejt1 = emissions[k,j,2:nObs]  # emissions for state j at t+1
                    bkwt1 = bkw[2:nObs,k,j] # backwards prob at t+1 for state j
                    #bkwt = bkw[1:(nObs - 1),k,i] # backwards prob at t for state i
                    fwdt = fwd[1:(nObs - 1),k,i] # forwards prob at t for state i

                    
                    #ot = ( postt * aij * bjt1 * bkwt1 ) / bkwt
                    ot = ( fwdt * aij * ejt1 * bkwt1 ) * distanceScale
                    #if(i!=j) os = ot/distanceScale  # scale back by distance between SNPs (except for diagonals, which are 1 - the rest)
                    #if(i==j) os = (ot - 1 + distanceScale)/distanceScale  # scale back by distance between SNPs (except for diagonals, which are different)                   
                    #return(os)
                    return(ot)
                },simplify="array")
                
            },simplify="array")
            
            o = apply(mat,c(2,3),sum) # sum over snps
            o = t(o)
            
            # note that estimates for diagonal values are not correct.
            return(o)
            
        },simplify="array")

    }

    # compute viterbi path as required
    viterbi = NULL
    if( computeViterbi ){
        
        viterbi = getViterbiPath(emissions,transitions)
        
    }
    
    print('return hmm')
    print(Sys.time())  

    return(list(posterior,transitionsPosterior,viterbi))
}



# path finding with emissions/transitions (currently a bit slow...)
getViterbiPath <- function(emissions,transitions){

    # set up arrays for viterbi
    nObs = dim(emissions)[3]
    nStates = dim(emissions)[2]   
    nInds = dim(emissions)[1]
    
    vertStates = array(dim = c(nObs,nInds,nStates)) # for saving viterbi values at each snp
    vertPath = array(dim = c(nObs,nInds))
    
    states = 1:nStates
    transitions_log = log(transitions)
    emissions_log = log(emissions)
    
    print("computing viterbi recursion...")

    for(i in 1:nObs){

        if(i%%1000==0) print(i)
        
        if(i==1)  {
            vertStates[i,,] = 0 # set starting state at arbitrarily 0
            maxValues = log(1/nStates) # probability of transition from starting to any of the the first states is 1/nstates.
        }

        if(i > 1) {
            
            toMaximise = sapply(states,function(l) vert_old + transitions_log[,l,i-1],simplify="array")
            maxValues = apply( toMaximise,c(2,3),function(x) max(x)) # what is the prob. value of the maximum state?
            maxStates = apply( toMaximise,c(2,3),function(x) which(x == max(x))) # which state maximises this prob?
            vertStates[i,,] = maxStates  # store maximum states for each transition for later.
            
        }
        
        vert_new = emissions_log[,,i] + maxValues
        vert_old = t(vert_new)
    }
    
    endState = apply( vert_old,2,function(x) which(x == max(x) ) )
    
    # back-tracking step
    print("back-tracing viterbi...")
    
    for(i in nObs:1){
        if(i == nObs) newState = endState     # set last postion to be the end state.        
        if(i < nObs) newState = sapply(1:nInds,function(ind) vertStates[i+1,ind,oldState[ind]])        
        vertPath[i,] = newState
        oldState = newState
    }

    return(vertPath)
}


# function for viterbi path finding from log2ratios (only use if not using myForwardsBackwards)
myViterbi <- function(log2Ratios,trans,scale_trans=NULL,HMMPosteriors){

    
    print('start hmm - viterbi')
    print(Sys.time())
                                        # trans = transition matrix (nStates x nStates)
                                        # log2Ratios = array:  nSNPs x nSamples
                                        # HMMPosteriors = array: nSNPs x nStates x 2
                                        # NOTE: The SNPs in Log2Ratios and HMMPosteriors must be in the same order!    
                                        # scale_a = an optional vector of length nObservations - 1 to reflect distance between SNPs
                                        # Values are in [0,1], where 1 is keep the baseline transition matrix
    
    
    # get transitions and emissions
    TE = getTransitionsAndEmissions(log2Ratios,trans,scale_trans,HMMPosteriors)
    
    transitions = TE$transitions
    emissions = TE$emissions        

    # get viterbi path
    vertPath = getViterbiPath(emissions,transitions)

    return(vertPath)
    
}



############
# likelihood equations
# AA and AB etc are functions to get likelihood given a gaussian. These must be prespecified in the R session

getLikelihoods <- function(Post,x){
        calls = dnorm(x, mean = Post[1], sd = Post[2], log = FALSE)
}


########################
# version of the hmm with multiple individuals at a time

# here, x is an nInds x 2 matrix of A and B intensities
getEmissions <- function(x,statePosteriors,ocean=0.0001){
  #x = the observation (log2R)
  #statePosteriors = estimates of mean and sds for this position.
  #emit = nSamples x nStates - emission probabilities for each state
    
    emit = sapply(1:dim(statePosteriors)[1],function(state){
        lhoods = getLikelihoods(statePosteriors[state,],x)
        l = (1-ocean) * lhoods + ocean
    })
    
    return(emit)
}
  

    
# get mean and variances for each SNP, at each copy number state used.

getSNPposteriors <- function(l2r,stateProbabilities,noise=FALSE,ocean=0.0001){

    # l2r = nSNPs x nSamples
    # stateProbabilities = nSNPS x nSamples x nStates
    # HMMPosteriors = nSNPs x nStates x 2 (the 2 comes from mean and variance)
    # NOTE: snp order in l2r and stateProbabilities must be the same!!!
    
    nSnps = dim(l2r)[1]
    nSamples = dim(l2r)[2]
    nStates = dim(stateProbabilities)[3]
    HMMPosteriors = array(dim=c(nSnps,nStates,2))

    print('computing means for each state')
    means = apply(stateProbabilities,3,function(stateProbs){
        # iterate over states (matrices nSNPs x nSamples)
        stateProbsFrac = stateProbs/rowSums(stateProbs)
        stateProbsFrac[is.na(stateProbsFrac)] = 0              
        m = rowSums(stateProbsFrac*l2r)
        return(m)
    })
    
    print('computing sds for each state')
    sds = apply(stateProbabilities,3,function(stateProbs){
                                        # iterate over states (matrices nSNPs x nSamples)
        s = sapply(1:nSnps,function(x) sqrt(wtd.var(l2r[x,],w=stateProbs[x,])))
        s[is.na(s)] = mean(s)  # impute sds where there is no information
        return(s)
    })

    # add a little gaussian noise in initialisation to avoid zeros in means (technically this value should be the background intensity... could estimate this )
    if( noise ) {
        print('adding some gaussian noise to means')
        means = means + rnorm(length(means),0,ocean)
    }
    
    HMMPosteriors[,,1] = means  # note SNPs are in the same order as in l2r
    HMMPosteriors[,,2] = sds
    
    return(HMMPosteriors)
}
    




# Function to get weighted variances. This is from package Hmisc    
wtd.var <- function (x, weights = NULL, normwt = TRUE, na.rm = FALSE, method = c("unbiased", 
    "ML")) 
{
    method <- match.arg(method)
    if (!length(weights)) {
        if (na.rm) 
            x <- x[!is.na(x)]
        return(var(x))
    }
    if (na.rm) {
        s <- !is.na(x + weights)
        x <- x[s]
        weights <- weights[s]
    }
    if (normwt) 
        weights <- weights * length(x)/sum(weights)
    if (method == "ML") 
        return(as.numeric(stats::cov.wt(cbind(x), weights, method = "ML")$cov))
    sw <- sum(weights)
    xbar <- sum(weights * x)/sw
    sum(weights * ((x - xbar)^2))/(sw - (if (normwt) 
        sum(weights^2)/sw
    else 1))
}


runHMM.em.OLD <- function(log2Ratios,initialStateProbabilities,maxIterations=3,trans=NULL,outIterations=2,converge.check=TRUE,tvdTol=1e-4,Expectation=NULL,ExpectationOrderOfSNPs=NULL,saveParameters=FALSE){

    # log2Ratios = nSNPs x nSamples
    # initialStateProbabilities = nSNPs x nSamples x nStates
    # NOTE: The order of SNPs in log2Ratios and initialStateProbabilities must be the same!
    
    print( paste0("HMM begun at ",date()) )
    postSNPs = rownames(log2Ratios)
    snpPositions = ps2snp$Position[match(postSNPs,ps2snp$ProbeSetID)] # positions for order in HMMPosteriors

    scale_trans = exp(-diff(sort(snpPositions))/10000)  # exponential decay in 10 KBs differences=    
    orderOfHMMSNPs = postSNPs[order(snpPositions)]   
    orderOfHMMSamples = colnames(log2Ratios)

    # order SNPs here by position.
    # order log2ratios by snp position
    orderedL2R = log2Ratios[order(snpPositions),]

    if(is.null(Expectation)) stateProbabilities = initialStateProbabilities[order(snpPositions),,]    # Note that if using existing Expectation values, then this must be in the same order as snps in log2Ratios.
        
    # check that parameter estimates are in the same order as orderedL2R
    if(!is.null(ExpectationOrderOfSNPs)) {
        Expectation = Expectation[match(rownames(orderedL2R),ExpectationOrderOfSNPs),,]
    }    

    HMM = list()
    convTvd = c()
    convCert = c()
    tvds = c()
    lastIteration = FALSE
    iteration = 0
#    for( iteration in 1:iterations){
    while( !lastIteration ){
        if(!is.null(Expectation)) lastIteration = TRUE  # just do one iteration if model parameters (Expectation) is provided.
            
        iteration = iteration + 1
        print( paste0("Iteration ",iteration) )
                                        # compute parameters (means, variances) for each of the snps and states
        if(is.null(Expectation)) Expectation = getSNPposteriors(l2r=orderedL2R,stateProbabilities)
                                        # compute posterior state probabilities
        Maximisation = myForwardsBackwards(orderedL2R,trans,trans_0=NULL,scale_trans=scale_trans,HMMPosteriors=Expectation)

        hmmOut = array(unlist(Maximisation),dim=c(nrow(Maximisation[[1]]),ncol(Maximisation[[1]]),length(Maximisation)))  # this is an inefficient step, just to change the dimension
        stateProbabilities = aperm(hmmOut, c(1,3,2))

        # compute some convergence statistics
        print('computing convergence stats')
        certs = compute.cert(stateProbabilities)
        convCert = c(convCert,mean(certs))
        
        if( iteration > 1 ){
            hmm.new = stateProbabilities
            tvds = compute.tvd(list(hmm.old,hmm.new))
            meanTvd = mean(tvds)
            convTvd = c(convTvd,meanTvd)
            if(( converge.check ) & (meanTvd < tvdTol) ) {
                lastIteration=TRUE
                print( paste0("convergence criteria (mean tvd < ",tvdTol,") met at iteration ",iteration) )
            }
        }
        hmm.old = stateProbabilities

                # always keep the last outIterations iterations until lastIteration==TRUE
        if( iteration <= outIterations ) HMM[[iteration]] = stateProbabilities
        if( iteration > outIterations ) {
            HMM = HMM[-1] # remove oldest version
            HMM[[outIterations]] = stateProbabilities # add in the new one
        }

        if( ( iteration == maxIterations ) & ( !lastIteration ) ) {
            print( paste0("convergence criteria not yet met at ",maxIterations," iterations. Stopping now anyway."))
            lastIteration=TRUE
        }        
    }

    # compile convergence statistics
    convStats = cbind(c(NA,convTvd),convCert)
    conv = list(tvds,certs)
    # Save parameter estimates that were used in the final iteration
    if(saveParameters) parameters = Expectation else parameters=NULL

    return(list("HMM"=HMM,"orderOfHMMSNPs"=orderOfHMMSNPs,"orderOfHMMSamples"=orderOfHMMSamples,"conv"=conv,"convStats"=convStats,"parameters"=parameters))
}    



runHMM.em <- function(log2Ratios,initialStateProbabilities,maxIterations=3,trans=NULL,outIterations=2,converge.check=TRUE,tvdTol=1e-4,Expectation=NULL,ExpectationOrderOfSNPs=NULL,saveParameters=TRUE,computeTransitions=TRUE){

    # log2Ratios = nSNPs x nSamples
    # initialStateProbabilities = nSNPs x nSamples x nStates
    # NOTE: The order of SNPs in log2Ratios and initialStateProbabilities must be the same!
    
    print( paste0("HMM begun at ",date()) )
    postSNPs = rownames(log2Ratios)
    snpPositions = ps2snp$Position[match(postSNPs,ps2snp$ProbeSetID)] # positions for order in HMMPosteriors

    scale_trans = exp(-diff(snpPositions[order(snpPositions)])/10000)  # exponential decay in 10 KBs differences=    
    orderOfHMMSNPs = postSNPs[order(snpPositions)]   
    orderOfHMMSamples = colnames(log2Ratios)

    # order SNPs here by position.
    # order log2ratios by snp position
    orderedL2R = log2Ratios[order(snpPositions),]
    
    doIterations = FALSE # by default, don't do iterations
    
    if(is.null(Expectation)) {
        doIterations = TRUE
        stateProbabilities = initialStateProbabilities[order(snpPositions),,]    # Note that if using existing Expectation values, then this must be in the same order as snps in log2Ratios.
    }
    
    # check that parameter estimates are in the same order as orderedL2R
    if(!is.null(ExpectationOrderOfSNPs)) {
        Expectation = Expectation[match(rownames(orderedL2R),ExpectationOrderOfSNPs),,]
    }    

    HMM = list()
    convTvd = c()
    convCert = c()
    tvds = c()
    lastIteration = FALSE
    iteration = 0
    initialExpectation = NULL
#    for( iteration in 1:iterations){
    while( !lastIteration ){
        if(!doIterations) lastIteration = TRUE  # just do one iteration if model parameters (Expectation) is provided.
            
        iteration = iteration + 1
        print( paste0("Iteration ",iteration) )
        
                                        # compute parameters (means, variances, transitions) for each of the snps and states
        if(doIterations) {
            # update emission parameters (actually this is maximisation step in EM!!)
            Expectation = getSNPposteriors(l2r=orderedL2R,stateProbabilities)
            
            if(iteration == 1) {
                print('Setting initial parameters all the same for each snp...')
            # If this is the first iteration,then set Expectations exactly the same for all SNPs + some random noise to the means. This avoids problems with bias in initial values. This was used only in version v3.                
                meanExp = apply(Expectation,c(2,3),mean) # take the mean of means and variances across all SNPs
                newExp = replicate(dim(Expectation)[1],meanExp)
                newExp[,1,] = newExp[,1,] + rnorm(dim(Expectation)[1]*dim(Expectation)[2],0,0.05) # random noise dim nStates x nSNPs
                Expectation = aperm(newExp, c(3,1,2))
                initialExpectation = Expectation
            }            
        }
        
       # update transition matrix (only possible after initial round of forwards/backwards)
        if( ( iteration > 1 ) & ( computeTransitions ) ) {

            print("updating transition probabilities...")
                                        # sum over transition matrices across all individuals.
            transPosteriors = Maximisation[[2]]
            posterior = Maximisation[[1]]
            newTrans = apply(transPosteriors,c(1,2),sum)/apply(posterior[1:(dim(posterior)[1]-1),,],2,sum)
            
                                        # diagonals make rows sum to 1
            for( i in 1:nrow(newTrans) ) newTrans[i,i] =  1 - sum(newTrans[i,-i])                
            trans = newTrans
            print(trans)
        }
        
                                        # compute posterior state probabilities (Actually this is the 'expectation' step in EM!!)
       # if(lastIteration) computeViterbi=TRUE else computeViterbi=FALSE # compute the viterbi path for the last iteration
        computeViterbi=FALSE
        Maximisation = myForwardsBackwards(orderedL2R,trans,trans_0=NULL,scale_trans=scale_trans,HMMPosteriors=Expectation,computeTransitions=computeTransitions,computeViterbi=computeViterbi)        

        
         # get new state probabilities (this is an inefficient step, just to change the dimensions!)
        stateProbabilities = aperm(Maximisation[[1]], c(1,3,2))

        # compute some convergence statistics
        print('computing convergence stats')
        certs = compute.cert(stateProbabilities)
        convCert = c(convCert,mean(certs))
        
        if( iteration > 1 ){
            hmm.new = stateProbabilities
            tvds = compute.tvd(list(hmm.old,hmm.new))
            meanTvd = mean(tvds)
            convTvd = c(convTvd,meanTvd)
            if(( converge.check ) & (meanTvd < tvdTol) ) {
                lastIteration=TRUE
                print( paste0("convergence criteria (mean tvd < ",tvdTol,") met at iteration ",iteration) )
            }
        }
        hmm.old = stateProbabilities

                # always keep the last outIterations iterations until lastIteration==TRUE
        if( iteration <= outIterations ) HMM[[iteration]] = stateProbabilities
        if( iteration > outIterations ) {
            HMM = HMM[-1] # remove oldest version
            HMM[[outIterations]] = stateProbabilities # add in the new one
        }

        if( ( iteration == maxIterations ) & ( !lastIteration ) ) {
            print( paste0("convergence criteria not yet met at ",maxIterations," iterations. Stopping now anyway."))
            lastIteration=TRUE
        }        
    }

    # compile convergence statistics
    convStats = cbind(c(NA,convTvd),convCert)
    conv = list(tvds,certs)
    
    # Save parameter estimates that were used in the final iteration
    if(saveParameters) parameters = list(Expectation,trans,initialExpectation) else parameters=NULL

    # extract viterbi path for last HMM version
    viterbi = Maximisation[[3]]
    
    return(list("HMM"=HMM,"viterbi"=viterbi,"orderOfHMMSNPs"=orderOfHMMSNPs,"orderOfHMMSamples"=orderOfHMMSamples,"conv"=conv,"convStats"=convStats,"parameters"=parameters))
    
}    




# convergence statistics
compute.tvd <- function(HMMpair){
    # HMM is list with sets of state probabilities.
    tvds = sapply(1:dim(HMMpair[[1]])[2],function(i){
        X = HMMpair[[1]][,i,]
        Y = HMMpair[[2]][,i,]  # this is the newest one
        tvd = 0.5*rowSums(abs(X-Y))
    })
    return(tvds)
}
compute.cert <- function(hmm){
    n = dim(hmm)[3]
    certs = sapply(1:dim(hmm)[2],function(i){
                                        # this computes how close each state probabilities are to c(0,1,0) etc.
        latest = hmm[,i,]
        cert = rowSums((latest-(1/n))^2)  # this is equivalent to length(a) * var(a)
    })
    certs2  = (n/(n-1))*certs
    return(certs2)
}


###############
# Read hdf5 version of hmm in batches
read.hmm.hdf5 <- function(h5,inds,snps=NULL,batch=NULL,otherInfo=NULL){

    if( (is.null(otherInfo)) & (is.null(batch)) ) {
        print("Can't extract. Must specify otherInfo or batch")
        return(NULL)
    }
                                        # snps must be in probesetID format
    if( (is.null(batch)) & (!is.null(otherInfo)) ){

        if(is.null(inds)) {
            print("Are you sure you want to extract ALL the cnv data??")
            theseBatches = unique(otherInfo$Batch)
        } else {
            theseBatches = unique(otherInfo$Batch[otherInfo$PIID %in% inds])
        }
        
    }
    
    if(!is.null(batch)) theseBatches=batch
    if(is.null(snps)) print('extracting all snps in the hmm')

    print("Reading hdf5 file")    
    batchesInH5 = unique(sapply(h5ls(h5)$name,function(s) str_split(s,pattern="\\.")[[1]][1]))
    missingBatches = theseBatches[!theseBatches%in%batchesInH5]
    if(length(missingBatches) > 0 ) {
        print( paste0("WARNING: the following batches are not in the hdf5 file. ",paste(missingBatches,collapse=",")) )
        theseBatches = theseBatches[theseBatches%in%batchesInH5]
    }


    # get snp orders (same for every batch)
    s1 = h5read(h5,paste0(theseBatches[1],".orderOfHMMSNPs"))
    if(!is.null(snps)) indexSnp = which(s1%in%snps) else indexSnp=c(1:length(snpOrder))
    snpOrderSubset = s1[indexSnp]

    
    # get sample orders
    dataBatchesSampleOrder = sapply(theseBatches,function(batch){
        
        sampleOrder = h5read(h5,paste0(batch,".orderOfHMMSamples"))

        if(!is.null(inds)) indexSamples = which(sampleOrder%in%inds) else indexSamples=c(1:length(sampleOrder))

        if((!is.null(inds)) & ( sum(inds%in%sampleOrder) == 0 ) ) {
            print("No relevant samples in this batch.")
            return(NULL)
        } else {
            return(sampleOrder[indexSamples])
        }
    },simplify=FALSE,USE.NAMES=FALSE)

    
    # get hmm results
    dataBatches = sapply(theseBatches,function(batch){
        print(batch)
        snpOrder = h5read(h5,paste0(batch,".orderOfHMMSNPs"))
        sampleOrder = h5read(h5,paste0(batch,".orderOfHMMSamples"))
        
        if(!is.null(snps)) indexSnp = which(snpOrder%in%snps) else indexSnp=c(1:length(snpOrder))
        if(!is.null(inds)) indexSamples = which(sampleOrder%in%inds) else indexSamples=c(1:length(sampleOrder))

        if((!is.null(inds)) & ( sum(inds%in%sampleOrder) == 0 ) ) {
            print("No relevant samples in this batch.")
            dataSubset=NULL
        } else {
            dataSubset = h5read(h5,paste0(batch,".HMM"),index=list(indexSnp,indexSamples,NULL))
        }
        return(dataSubset)
    },simplify=FALSE)

    parameters = h5read(h5,"parameters")
    #parameters = parameters[indexSnp,,]
    parameters = parameters[[1]][indexSnp,,]
    
    # create list as plotting expects
    hmm = list()
    hmm[[1]] = abind(dataBatches,along = 2)
    orderOfHMMSNPs = snpOrderSubset
    orderOfHMMSamples = abind(dataBatchesSampleOrder,along=1)
    
    outhmm  = list("HMM" = hmm,"orderOfHMMSNPs"=orderOfHMMSNPs,"orderOfHMMSamples"=orderOfHMMSamples,"parameters"=parameters)
    
    return(outhmm)    
}




###############
# functions for plotting hmm
    
transformHMMtoCalls <- function(hmmToPlot,samplesToPlot=NULL,snpsToPlot=NULL,threshold){
# this function makes hard calls of copy number from hmm probabilities    
# hmmToPlot is output from runHMM.em. pull out the last iteration.
    
    data = hmmToPlot$HMM[[length(hmmToPlot$HMM)]]  # last in HMM list is the last iteration
    snps = outRownames = hmmToPlot$orderOfHMMSNPs
    sampleIDs = outColnames = hmmToPlot$orderOfHMMSamples
    nStates = dim(data)[3]
    
    if( threshold == "viterbi" ) data = hmmToPlot$viterbi
        
   # subsets for snps and samples and create matrix
       
    if(!is.null(snpsToPlot)) {
        if( length(snpsToPlot)==1 ) snpsToPlot = c(snpsToPlot,snpsToPlot)
        snpSubset = snps%in%snpsToPlot       
        if( threshold == "viterbi" ) data = data[snpSubset,] else data = data[snpSubset,,]
        outRownames = snps[snpSubset]
    }
    
    if(!is.null(samplesToPlot)) {
        
        if( length(samplesToPlot)==1 ) samplesToPlot = c(samplesToPlot,sampleIDs[sampleIDs!=samplesToPlot][1])

        sampleSubset = sampleIDs%in%samplesToPlot
        if( threshold == "viterbi" ) data = data[,sampleSubset] else data = data[,sampleSubset,]
        outColnames = sampleIDs[sampleSubset]
    }

    if( threshold != "viterbi" ) {
                                        # record states above threshold
        a = data >= threshold
        
        forPlot = a[,,1]    
        for(s in 2:nStates){        
            forPlot = forPlot + s*a[,,s]
        }        
                                        #   forPlot = apply(a,2,function(b) b%*%c(1:nStates))   # this is slower
    }

    if( threshold == "viterbi" ) forPlot = data
        
    colnames(forPlot) = outColnames
    rownames(forPlot) = outRownames
    
    return(forPlot)
    
}


hclust.order <- function(W){

    Wna = W; Wna[Wna==0] = NA
    d = dist(t(Wna)) # euclidian distance, ignoring uncertain values coded as NAs. takes a few seconds on 1500 x 1500
    j = hclust(d)
    hclustOrderSamples = colnames(W)[j$order]
    return(hclustOrderSamples)
}

colors = c("lightgray","blue","red","green") # for 3 state HMM
names(colors) = as.character(0:3)

getXY <-  function(W,snpPositions,cols = colors,scale=FALSE,colLims=NULL,threshold=NULL,colorSet=c("red","blue","green"),nCols=10){
    # colors should be length of unique(c(W)), if W is integer
# W is matrix nSNPs x nSamples, with snp and sample names attached
    values = sort( unique(c(0,c(W))) )  # include 0 by default
    if( length(values) > 20 ) {
        scale=TRUE
        print("too many unique values. using scaled colors.")
    } else {
        if( length(cols) < length(values) ){
            print("Not enough colors specified. Making them for you.")
            choices = colors()[!grepl("gray|white|snow",colors())]
            cols = sample(choices,length(values),replace=FALSE)
            names(cols)=values
        }
    }
    
    y = unlist(lapply(1:ncol(W),function(j) rep(j,dim(W)[1])))
    p = unlist(lapply(1:ncol(W), function(i) W[,i]))
    x=rep(snpPositions,dim(W)[2])

    if(scale==F){  # values are hard calls
        cols3 = rep(NA,length(p))
        cols3 = cols[as.character(p)]
    } else {
        
        if(is.null(colLims)) colLims = c(-max(abs(W),na.rm=T),max(abs(W),na.rm=T))
        
        cols3 = rep("gray",length(p))
        cols3 = GetColourScale2(p,colourSet=colorSet,nBreaks=nCols,fixedLims=c(-threshold,threshold))
        cols3[is.na(cols3)] = "gray"                
    }
    return(list("x"=x,"y"=y,"cols"=cols3,"colLabels"=cols))
}


plotHMM <- function(hmmToPlot,samplesToPlot,snpsToPlot=NULL,filename,threshold=0.9,title="",plotRecomb=TRUE,extraPlotting=FALSE,Width=2000,Height=1000,recombRates=NULL,sortInds=FALSE,updateLegend=NULL,showPar=TRUE,pngRes=150,addGenes=FALSE,spaceForGenes=0.2,...){

    W = transformHMMtoCalls(hmmToPlot,samplesToPlot,snpsToPlot,threshold)
    if(sortInds) {
        Wna = W; Wna[Wna==0] = NA
        # by default - order by fractions of each 'color'
#        cats = unique(c(Wna))        
#        k = sapply(sort(cats[!is.na(cats)]),function(cat) k = colSums(Wna==cat,na.rm=TRUE)/colSums(!is.na(Wna)) )
#        Order = do.call(order, as.data.frame(k))
        Order = order(colMeans(Wna,na.rm=TRUE))
    } else {
        if( sum(!samplesToPlot%in%colnames(W)) > 0 ) print('Not all samples requested are actually in the input data. Removing these from requested list.')
        samplesToPlot = samplesToPlot[samplesToPlot%in%colnames(W)]
        Order = match(samplesToPlot,colnames(W))  # order by input: samplesToPlot        
    }

    snpPositions = ps2snp$Position[match(rownames(W),ps2snp$ProbeSetID)]
    plottingValues = getXY(W[,Order],snpPositions)
    legendLabels = names(plottingValues$colLabels)
    if(!is.null(updateLegend)) legendLabels[match(names(updateLegend),legendLabels)] = unlist(updateLegend)
        
    offSet= max(plottingValues$y)/20
                                        #par(bg="lightgray")
    print("Creating plot")
    
    png(paste0(filename,'-threshold',threshold,'.png'),width=Width,height=Height,res=pngRes)
        
    if(addGenes) {
        layout(mat=c(1,2),heights=c( 1-spaceForGenes,spaceForGenes ) )
        extraPlotting=TRUE
    }
    
    plot(plottingValues$x,plottingValues$y,col=plottingValues$cols,pch=15,xlab="position on X (MB)",ylab="individuals",ylim=c(-offSet,max(plottingValues$y)),yaxt="n",xaxt="n",main=paste0(title,"\nHMM posterior >= ",threshold),...)
    axis(2,at=seq(0,max(plottingValues$y),by=100))
    axis(1,at = axTicks(1),labels=axTicks(1)/100000)
    
    legend("topleft",legend=legendLabels,col=plottingValues$colLabels,pch=15,bty="n")
    
    if(plotRecomb==T){
        lines(recombRates$position,recombRates$COMBINED_rate.cM.Mb./(max(recombRates$COMBINED_rate.cM.Mb.)/offSet) - offSet,col="black",xlim=c(min(plottingValues$x,na.rm=T),max(plottingValues$x,na.rm=T)))
    }

    if(showPar){
        abline(v=par1,lty=3)
        text(x=par1,y=par()$usr[3],labels="PAR1",pos=1,lwd=3,xpd=NA)
        abline(v=par1,lty=3)
        text(x=par2,y=par()$usr[3],labels="PAR2",pos=1,lwd=3,xpd=NA)
    }
    
    if(!extraPlotting) dev.off() else print("Don't forget to run dev.off() eventually!")
    
}

# Function to plot parameters
plotParameters <- function(hmmToPlot,colors,snpsToPlot,add=FALSE,ylimits=NULL,xlimits=NULL,...){

    snpIndex = match(snpsToPlot,hmmToPlot$orderOfHMMSNPs)
    snpPositions = ps2snp$Position[match(snpsToPlot,ps2snp$ProbeSetID)]
    modelParameters = hmmToPlot$parameters#[[1]]
        
    means = modelParameters[snpIndex,,1]
    sds = modelParameters[snpIndex,,2]
    #quant = qnorm(0.025,0,1)
    quant = 1 # just one standard deviation
    quant975 = means + quant*sds
    quant025 = means - quant*sds
    maxQuant = max( c(max(abs(quant975)),max(abs(quant025))) )
    
    x = snpPositions
    if(is.null(xlimits)) xlimits = c(min(x),max(x))
    if(is.null(ylimits)) ylimits = c(-maxQuant,maxQuant)

    if(!add) plot(NULL,xlim=xlimits,ylim=ylimits,ylab="log2Ratio",xlab="Position",...)

    for( state in 1:ncol(means)){        
        color = colors[as.character(state)]
        lines(x,means[,state],col=color,lwd=0.3)
        lines(x,quant975[,state],lty=3,col=color,lwd=0.3)
        lines(x,quant025[,state],lty=3,col=color,lwd=0.3)        

    }    

}

plotParametersHist <- function(hmmToPlot,modelParameters=NULL,colors=NULL,snpsToPlot=NULL,ylimits=NULL,xlimits=NULL,...){
    
    if( is.null(snpsToPlot) & (!is.null(hmmToPlot)) ) snpsToPlot = hmmToPlot$orderOfHMMSNPs

    if( !is.null(hmmToPlot) ) {
        snpIndex = hmmToPlot$orderOfHMMSNPs%in%snpsToPlot
        modelParameters = hmmToPlot$parameters[[1]]
    } else {
        snpIndex = rep(TRUE,dim(modelParameters)[1])
    }
    
    means = modelParameters[snpIndex,,1]
    sds = modelParameters[snpIndex,,2]
    variances = sds^2
    nSnps = length(snpIndex)
    
    par(mfrow=c(ncol(means),2),mar=c(2,4,1,1))    
    xlimMean = c(min(means),max(means))
    xlimVar = c(min(variances),max(variances))
    
    for( state in 1:ncol(means)){
        if( state == 1 ){
            t1="Means"
            t2="Variances"
        } else {
            t1 = t2 = NA
        }
        if( is.null(colors) ) color = "black" else color = colors[as.character(state)]
        
        hist(means[,state],xlim=xlimMean,main=t1,ylab=state,breaks=nSnps/100,border=NA,col=color)
        hist(variances[,state],xlim=xlimVar,main=t2,ylab=NA,breaks=nSnps/100,border=NA,col=color)
    }    

}


# Function to plot the l2r of one individual, along with the mean/variance estimates from a specific hmm run.
plotL2RIndividual <- function(log2Ratios,hmmToPlot,individual,snpsToPlot=NULL,filename,threshold=0.9,title="",plotRecomb=TRUE,extraPlotting=FALSE,Width=2000,Height=1000,recombRates=NULL,updateLegend=NULL,showPar=TRUE,xlimits=NULL,ylimits=NULL,modelParameters=NULL,uncertainColour="gray",...){

    #hmmToPlot = dataHMM86
    #individual = samplesToPlotKEEP[1700]
    #snpsToPlot = theseSnps
    #threshold=0.98
    #log2Ratios = read.cnv.hdf5(inds=samplesToPlotKEEP,snps=snps,type="log2ratio")
    path = transformHMMtoCalls(hmmToPlot,individual,snpsToPlot,threshold)
    path1 = path[,individual]
    snps = names(path1)
    
    l2r = log2Ratios[snps,individual]
    snpPositions = ps2snp$Position[match(snps,ps2snp$ProbeSetID)]
    
    # get values to plot
    plottingValues = getXY(path,snpPositions)
    x = plottingValues$x[-(1:(length(plottingValues$x)/2))]    
    plottingValues$cols = plottingValues$cols[-(1:(length(plottingValues$x)/2))]
    
    if( sum(as.numeric(names(plottingValues$cols))!=path1) > 0 ) print("WRONG!!!")
    
    legendLabels = names(colors)
    if(!is.null(updateLegend)) legendLabels[match(names(updateLegend),legendLabels)] = unlist(updateLegend)

    if( is.null(modelParameters) ) modelParameters = hmmToPlot$parameters#[[1]]
    if( is.null(modelParameters) ){
        print("ERROR: no model parameters for emissions are provided.")
    }
    
    png(paste0(filename,'-threshold',threshold,'-',individual,'.png'),width=Width,height=Height,res=150)

    # plot l2r raw data
    x = snpPositions
    y = l2r

    if(is.null(xlimits)) xlimits = c(min(x),max(x))
    if(is.null(ylimits)) ylimits = c(-max(abs(y)),max(abs(y)))
    

    # plot means and variances used for each snp
    plotParameters(hmmToPlot,colors=colors,xlimits=xlimits,snpsToPlot,add=FALSE,main=paste0(title,"\nHMM posterior >= ",threshold))

    points(x,y,col=plottingValues$cols,pch=16,cex=0.5)

                                        # highligh path for this individual
    states = as.character(c(0,1,2,3))
    paths = sapply(states,function(state) path1==state )
    nonZeroPaths = paths[,colnames(paths)!="0",drop=FALSE]
    snpIndex = match(snpsToPlot,hmmToPlot$orderOfHMMSNPs)
    means = modelParameters[snpIndex,,1]
    
    y2 = rowSums(means[,as.numeric(colnames(nonZeroPaths))]*nonZeroPaths)
#    lineColours = plottingValues$cols
#    lineColours[y2==0] = "transparent"
    y2[y2==0] = NA
    
    lines(x,y2)

    legend("topleft",legend=legendLabels,col=colors,pch=15,bty="n")
    
    if(plotRecomb==TRUE){
        offSet = (par()$usr[4] - par()$usr[3])/10

        lines(recombRates$position,par()$usr[3]+recombRates$COMBINED_rate.cM.Mb./(max(recombRates$COMBINED_rate.cM.Mb.)/offSet),col="black",xlim=c(min(x,na.rm=T),max(x,na.rm=T)))
    }

    
    if(showPar){
        abline(v=par1,lty=3)
        text(x=par1,y=par()$usr[3],labels="PAR1",pos=1,lwd=3,xpd=NA)
        abline(v=par1,lty=3)
        text(x=par2,y=par()$usr[3],labels="PAR2",pos=1,lwd=3,xpd=NA)
    }
    
    if(!extraPlotting) dev.off() else print("Don't forget to run dev.off() eventually!")

}
