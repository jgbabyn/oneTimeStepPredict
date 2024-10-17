#' Create an index mapping between the parameter vector and time steps
#'
#' This function maps what time steps correspond to which parameter in the list
#'
#' @param pars this is the vector of parameters from the new object
#' @param data.vector.indicators the names of the DATA INDICATORS in the TMB template
#' @param observation.names the names corresponding to the observations
#' @param time.step.indicators the list of vectors giving the time steps to match up
#' @param obs_support the support of the observations
#' @param c_or_d whether or not the observation is continuous or discrete 
create_t_index_key <- function(pars,data.vector.indicators,observation.names,time.step.indicators,obs_support,c_or_d,cdf=FALSE){
    tdf = data.frame(name=names(pars))
    tdf$oind = 1:nrow(tdf)
    tdf$ts = NA
    tdf$c_or_d = NA
    tdf$obs_LL = NA
    tdf$obs_UL = NA
    for(i in 1:length(data.vector.indicators)){
        tdf$ts[which(tdf$name == data.vector.indicators[i])] = time.step.indicators[[i]]
        tdf$ts[which(tdf$name == observation.names[i])] = time.step.indicators[[i]]
        tdf$c_or_d[which(tdf$name == observation.names[i])] = c_or_d[i]
        tdf$c_or_d[which(tdf$name == data.vector.indicators[i])] = c_or_d[i]
        tdf$obs_LL[which(tdf$name == observation.names[i])] = obs_support[[i]][1]
        tdf$obs_UL[which(tdf$name == observation.names[i])] = obs_support[[i]][2]
        tdf$obs_LL[which(tdf$name == data.vector.indicators[i])] = obs_support[[i]][1]
        tdf$obs_UL[which(tdf$name == data.vector.indicators[i])] = obs_support[[i]][2]
        if(cdf == TRUE){
            tdf$ts[which(tdf$name == paste0(data.vector.indicators[i],".lower_cdf"))] = time.step.indicators[[i]]
            tdf$ts[which(tdf$name == paste0(data.vector.indicators[i],".upper_cdf"))] = time.step.indicators[[i]]
        }
 
    }
    tdf
}


#' Get the ZCA whitening matrix
#'
#' @param Hess the hessian matrix of the obs
#' @param just_cor get just the correlation whitening matrix?
#'
ZCA_whiten <- function(Hess,just_cor = FALSE){
    if(just_cor == TRUE){
        sigma = solve(Hess)
        cory = cov2cor(sigma)
        invsigsq = expm::sqrtm(solve(cory))
    }else{
        invsigsq = expm::sqrtm(Hess)
    }
    invsigsq
}

#' Get the Cholesky Whitening Matrix
#'
#' @param Hess the hessian matrix of the observation
#' @param just_cor get just the correlation whitening matrix?
#' @param no_inv don't do the inverse?
#'
chol_whiten <- function(Hess,just_cor = FALSE, no_inv=FALSE){
    Sigma = solve(Hess)
    if(just_cor == FALSE){
        L = t(chol(Sigma))
        if(no_inv == FALSE){
            invL = solve(L)
        }else{
            invL = L
        }
        
    }else{
        cory = cov2cor(Sigma)
        if(no_inv == FALSE){
            invL = solve(t(chol(cory)))
        }else{
            invL = chol(cory)
        }
    }
    invL
}

#' Get the ZCA-Correlation Whitening Matrix
#'
#' @param Hess the hessian matrix of the observation
#' @param just_cor get just the correlation whitening matrix?
ZCACor_whiten <- function(Hess,just_cor = FALSE){
    Sigma = solve(Hess)
    CorrMI = solve(cov2cor(Sigma))
    PinvSq = expm::sqrtm(CorrMI)
    if(just_cor == FALSE){
        VinvSq = diag(1/sqrt(diag(Sigma)))
        W = PinvSq%*%VinvSq
    }else{
        W = PinvSq
    }
    W
}


#' Get the indices corresponding to a given time step
#'
#' @param ts the time step to get the indices corresponding to
#' @param vname the vector of names to match to
#' @param t_index the index mapping between parameter and time step (from create_t_index_key)
#' @param grthan return the values greater than the time step?
get_yinds <- function(ts,vname,t_index,grthan=FALSE){
    if(grthan == FALSE){
        t_index$oind[(t_index$name %in% vname) & t_index$ts == ts]
    }else{
        t_index$oind[(t_index$name %in% vname) & t_index$ts > ts]
    }
}

#' Get the indices corresponding to a certain observation
#'
#' @param ts the certain observation to get the indices corresponding to
#' @param vname the vector of names to match to
#' @param t_index the index mapping between parameter and time step (from create_t_index_key)
#' @param grthan return the values greater than the time step?
get_lin_inds <- function(ts,vname,t_index,grthan=FALSE){
    if(grthan == FALSE){
        t_index$oind[(t_index$name %in% vname) & t_index$lin_order == ts]
    }else{
        t_index$oind[(t_index$name %in% vname) & t_index$lin_order > ts]
    }
}

#' Function to generate the function to control what variables are used during a time step
#'
#' Creates a function control what observations are included (non-zero weighted) during a time step as
#' well as set the current time step variables to a given vector of values.
#'
#' @param pars the vector of parameters from the new object
#' @param t_index the index mapping between parameter and time step (from create_t_index_key)
#' @param data.vector.indicators the names of the DATA INDICATORS in the TMB template
#' @param observation.names the names corresponding to the observations
toggle_obs_gen <- function(pars,t_index,data.vector.indicators,observation.names){

    toggle_obs <- function(ts,p){
        wt_inds = get_yinds(ts,data.vector.indicators,t_index,TRUE)
        pars[wt_inds] = 0
        obs_inds = get_yinds(ts,observation.names,t_index,FALSE)
        pars[obs_inds] = p
        pars
    }
    toggle_obs
}


#' Function to generate the function to control what variables are used during a time step
#'
#' Creates a function control what observations are included (non-zero weighted) during a time step as
#' well as set the current time step variables to a given vector of values.
#'
#' @param pars the vector of parameters from the new object
#' @param t_index the index mapping between parameter and time step (from create_t_index_key)
#' @param data.vector.indicators the names of the DATA INDICATORS in the TMB template
#' @param observation.names the names corresponding to the observations
lin_toggle_obs_gen <- function(pars,t_index,data.vector.indicators,observation.names){

    toggle_obs <- function(ts,p){
        wt_inds = get_lin_inds(ts,data.vector.indicators,t_index,TRUE)
        pars[wt_inds] = 0
        obs_inds = get_lin_inds(ts,observation.names,t_index,FALSE)
        pars[obs_inds] = p
        pars
    }
    toggle_obs
}

#' Find the residuals for one time step assuming the conditional distribution follows a normal distribution
#'
#' @param ts the current time step
#' @param obj the newly created obj with each observation as a parameter
#' @param data.vector.indicators the names of the DATA INDICATORS in the TMB template
#' @param observation.names the names corresponding to the observations
#' @param w_method the whitening method to use
oneTimeStepGaussian <- function(ts,obj,t_index,data.vector.indicators,observation.names,w_method){

    toggle_obs <- toggle_obs_gen(obj$par,t_index,data.vector.indicators,observation.names)
    ts_ind = get_yinds(ts,observation.names,t_index,FALSE)
 
    fn <- function(y){
        obj$fn(toggle_obs(ts,y))
    }

    gr <- function(y){
        obj$gr(toggle_obs(ts,y))[ts_ind]
    }
        
    obs <- obj$par[ts_ind]

    opt = nlminb(obs,fn,gr)
    cent = obs - opt$par
    Hess <- optimHess(opt$par,fn,gr)
    ## if(cholesky == FALSE){
    ##     invsigsq = expm::sqrtm(Hess)
    ##     pred_residual = t(cent)%*%invsigsq
    ## }else{
    ##     Sigma = solve(Hess)
    ##     L = t(chol(Sigma))
    ##     pred_residual = as.vector(solve(L,cent))
    ## }
    W = switch(w_method,
               ZCA=ZCA_whiten(Hess,FALSE),
               ZCACor=ZCACor_whiten(Hess,FALSE),
               Cholesky=chol_whiten(Hess,FALSE))

    ##pred_residual = t(cent)%*%W
    pred_residual = W%*%cent
    list(residuals=pred_residual,predictions=opt$par,raw_resid=cent)
}


oneTimeStepGeneric <- function(ts,obj,t_index,data.vector.indicators,observation.names,w_method){
    print(ts)

    toggle_obs <- toggle_obs_gen(obj$par,t_index,data.vector.indicators,observation.names)
    ts_ind = get_yinds(ts,observation.names,t_index,FALSE)
    obs = obj$par[ts_ind]
    mobs = matrix(obs,ncol=1)

    t_range = as.matrix(t_index[t_index$oind %in% ts_ind,c("obs_LL","obs_UL")])
    range = t(t_range)
    
    fn <- function(y){
        obj$fn(toggle_obs(ts,y))
    }
    nll = fn(mobs)

    ##Takes a matrix and returns a matrix
    vfn <- function(y){
        apply(y,2,function(x){obj$fn(toggle_obs(ts,x))})
    }
    

    gr <- function(y){
        obj$gr(toggle_obs(ts,y))[ts_ind]
    }
    opt = nlminb(obs,fn,gr)
    Hess <- optimHess(opt$par,fn,gr)
    
    F_fn <- function(y){
       exp(-(fn(y)-nll))
    }

    ##Takes a matrix and returns a matrix
    vF_fn <- function(y){
        matrix(apply(y,2,function(x){exp(-(fn(x)-nll))}),ncol=ncol(y))
    }
    ###try(UDT <- cubature::hcubature(vF_fn,lowerLimit=range[1,],upperLimit=range[2,],vectorInterface=TRUE))

    UD = cubature::hcubature(vF_fn,lowerLimit=range[1,],upperLimit=range[2,],vectorInterface=TRUE)$integral
    UN = numeric(length(obs))
    for(i in 1:length(UN)){
        trange = range[2,]
        trange[i] = obs[i]
        UN[i] = cubature::hcubature(vF_fn,lowerLimit=range[1,],upperLimit=trange,vectorInterface=TRUE)$integral
    }
    Us = UN/UD
    ##Under the normal distribution but still correlated
    inNorm = qnorm(Us)
     W = switch(w_method,
               ZCA=ZCA_whiten(Hess,TRUE),
               ZCACor=ZCACor_whiten(Hess,TRUE),
               Cholesky=chol_whiten(Hess,TRUE))

    
    pred_residual = W%*%inNorm
    ##pred_residual = t(inNorm)%*%W
    pred_residual
}

oneTimeStepWorkaround <- function(ts,obj,t_index,data.vector.indicators,observation.names,w_method,t_index2){

    toggle_obs <- toggle_obs_gen(obj$par,t_index,data.vector.indicators,observation.names)
    lin_toggle_obs <- toggle_obs_gen(obj$par,t_index2,data.vector.indicators,observation.names)

    
    ts_ind = get_yinds(ts,observation.names,t_index,FALSE)
    obs = obj$par[ts_ind]

    t_range = as.matrix(t_index[t_index$oind %in% ts_ind,c("obs_LL","obs_UL")])
    range = t(t_range)
    
    integral_res = numeric(length(ts_ind))
    cdf_lower = numeric(length(ts_ind))
    cdf_upper = numeric(length(ts_ind))

    
    
    ##First find the oneStepAhead residuals so we can recorrelate them
    for(i in 1:length(ts_ind)){
        
        c_ind = t_index2[ts_ind[i],c("ts")]
        c_or_d = t_index2[ts_ind[i],c("c_or_d")]
        i_fn <- function(y){
            obj$fn(lin_toggle_obs(c_ind,y))
        }

        i_gr <- function(y){
            obj$gr(lin_toggle_obs(c_ind,y))[c_ind]
        }

        nll = i_fn(obs[i])
        i_F <- function(y){
            exp(-(i_fn(y)-nll))
        }

        i_vF = Vectorize(i_F)

        do_continuous <- function(){
            ##Integrate works better for 1D?
            UD = try({integrate(i_vF,range[1,i],range[2,i])$value})
            UN = try({integrate(i_vF,range[1,i],obs[i])$value})
            #UD = cubature::hcubature(i_F,lowerLimit = range[1,i],upperLimit = range[2,i])$integral 
                                        #UN = cubature::hcubature(i_F,lowerLimit = range[1,i],upperLimit = obs[i])$integral
            if(is(UN,"try-error")){
                UN = NaN
            }
            if(is(UD,"try-error")){
                UD = NaN
            }
            UN/UD
        }

        do_discrete <- function(){
            support = seq(range[1,i],range[2,i],by=1)
            vF = Vectorize(i_F)(support)
            lowF = sum(vF[support <= obs[i]])
            highF = sum(vF[support > obs[i]])
            cdfL = nll - log(lowF)
            cdfU = nll - log(highF)
            Fx = 1/(1+exp(log(highF)-log(lowF)))
            px = 1/(exp(-log(lowF))+exp(-log(highF)))
            Fx-runif(1)*px
        }

        if(c_or_d == "continuous"){
            integral_res[i] = do_continuous()
        }else{
            integral_res[i] = do_discrete()
        }
        

    }

    ##Next find the correlations of the predictions at a time step
    fn <- function(y){
        obj$fn(toggle_obs(ts,y))
    }

    gr <- function(y){
        obj$gr(toggle_obs(ts,y))[ts_ind]
    }
    
    obs <- obj$par[ts_ind]

    opt = nlminb(obs,fn,gr)
    cent = obs - opt$par
    Hess <- optimHess(opt$par,fn,gr)
    
    zorp = try({W = switch(w_method,
               ZCA=ZCA_whiten(Hess,FALSE),
               ZCACor=ZCACor_whiten(Hess,FALSE),
               Cholesky=chol_whiten(Hess,FALSE))
    WChol = chol_whiten(Hess,FALSE,TRUE)
    n_corr_res = qnorm(integral_res)
    n_uncorr_res = WChol%*%n_corr_res
               pred_residual =W%*%n_uncorr_res})
    if(is(zorp,"try-error")){
        pred_residual = rep(NaN,length(integral_res))
    }
    ret = data.frame(residual = pred_residual,F=integral_res,ts=ts,obs_name=names(obs))
    ret

   
    
}



#' Find the standardized residuals predicting one time step at a time
#'
#' @param obj the optimized TMB model object to find the residuals of
#' @param data.vector.indicators the names of the DATA INDICATORS in the TMB template
#' @param observation.names the names corresponding to the observations
#' @param time.step.indicators vector or list of vectors giving the time steps corresponding to the observations/data.vector.indicators
#' @param method Which method to use, currently only oneTimeStepGaussian or oneTimeStepGeneric is supported
#' @param c_or_d vector indicating whether or not observations are continuous or discrete, continuous assumed by default
#' @param obs_support the list of the lower and upper support of the observations vectors (optional, assumes -\infty to \infty)
#' @param w_method choice of method to use for whitening variables, either 'ZCACor' (default), 'ZCA' and 'Cholesky'
#'
#' @export
oneTimeStepPredict <- function(obj,observation.names=NULL,data.vector.indicators=NULL,time.step.indicators=NULL,method=c("oneTimeStepGaussian"),c_or_d=rep("continuous",length(observation.names)),obs_support=NULL,w_method=c('ZCACor')){

    if(missing(obj)){
        stop("Need a TMB model object to calculate the residuals of!")
    }
    if(missing(observation.names)){
        stop("Need names of observations in the TMB template!")
    }
    if(missing(data.vector.indicators)){
        stop("Need names corresponding to data vector indicators in the TMB template!")
    }
    if(missing(time.step.indicators)){
        stop("Need vector or list of vectors indicating the time steps corresponding to each vector of observations!")
    }

    if(is.vector(time.step.indicators) & !is.list(time.step.indicators)){
        tmp_ts_list = list(time.step.indicators)
        time.step.indicators=tmp_ts_list
    }

    args <- as.list(obj$env)[intersect(formalArgs(TMB::MakeADFun),ls(obj$env))]

    if(!all(observation.names %in% names(args$data))){
        stop("Observation names must be in template!")
    }

    if(is.null(obs_support)){
        obs_support = list()
        for(i in 1:length(observation.names)){
            obs_support[[i]] = c(-Inf,Inf)
        }
    }
    
    
    
    unique_ts = unique(do.call(c,time.step.indicators))

    ##Check if Hessian is usable internally
    intheCheck = try(obj$he(),silent=TRUE)
    if(is(intheCheck,"try-error")){
        intHe = FALSE
    }else{
        intHe = TRUE
    }
    
    if(length(obj$env$random) > 0){
        args$parameters = obj$env$parList(par=obj$env$last.par.best)
    }else{
        args$parameters = obj$env$parList(obj$env$last.par.best)
    }
    parm_names = names(args$parameters)
    random_names = unique(obj$env$.random)
    fixed_names = setdiff(parm_names,random_names)
    o_map = args$map
    f_map = lapply(args$parameters[fixed_names], function(x) as.factor(x*NA))
    o_map[fixed_names] = f_map

    args$parameters[observation.names] = args$data[observation.names]
    args$parameters[data.vector.indicators] = lapply(time.step.indicators,function(x){rep(1,length(x))})
    args$data[observation.names] = NULL
    args$map = o_map
    args$random = random_names
    args$silent = TRUE
    if(intHe == TRUE & length(random_names) > 0){
        args$intern = TRUE
    }
    n_obj = do.call(TMB::MakeADFun,args)

    t_index = create_t_index_key(n_obj$par,data.vector.indicators,observation.names,time.step.indicators,obs_support,c_or_d,FALSE)
    og_lens = lapply(time.step.indicators,function(x){
        ret = 1:length(x)
    })
    t_index2 = create_t_index_key(n_obj$par,data.vector.indicators,observation.names,og_lens,obs_support,c_or_d,FALSE)
    
    
    if(method == "oneTimeStepGaussian"){
        residuals = n_obj$par
        predictions = n_obj$par
        for(y in unique_ts){
            ts_ind = get_yinds(y,observation.names,t_index,FALSE)
            residuals[ts_ind] = oneTimeStepGaussian(y,n_obj,t_index,data.vector.indicators,observation.names,w_method)$residuals
            predictions[ts_ind] = oneTimeStepGaussian(y,n_obj,t_index,data.vector.indicators,observation.names,w_method)$predictions
        }
        t_index$residuals = residuals
        t_index$predictions = predictions
        ret = t_index[t_index$name %in% observation.names,c("name","ts","residuals","predictions")]
    }
    if(method == "oneTimeStepGeneric"){
        residuals = n_obj$par
        predictions = n_obj$par
        for(y in unique_ts){
            ts_ind = get_yinds(y,observation.names,t_index,FALSE)
            residuals[ts_ind] = oneTimeStepGeneric(y,n_obj,t_index,data.vector.indicators,observation.names,w_method)
        }
        t_index$residuals = residuals
        ret = t_index[t_index$name %in% observation.names,c("name","ts","residuals")]

    }
    if(method == "oneTimeStepAny"){
        rom = list()
        for(j in 1:length(unique_ts)){
            y = unique_ts[j]
            rom[[j]] = oneTimeStepWorkaround(y,n_obj,t_index,data.vector.indicators,observation.names,w_method,t_index2)
        }
        romC = do.call(rbind,rom)
        ret = romC
    }
        

    ret
}
