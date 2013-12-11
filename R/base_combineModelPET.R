# combine model with same dimension

.combineModelPET <- function( fit_list, mu_list, pi_list, pi0_list,
    gamma_list, BIC_vec, AIC_vec, maxComp ) { 
    
    fit_list_new <- mu_list_new <- pi_list_new <- pi0_list_new <- 
        gamma_list_new <- vector( "list", maxComp )        
    BIC_vec_new <- AIC_vec_new <- rep( NA, maxComp )    
    dimVec <- unlist(lapply( mu_list, length ))
    
    for ( i in 1:maxComp ) {
        model_i <- which( dimVec==i )
        
        if ( length(model_i) > 0 ) {
            # if more than one model have the same dim, choose one with smaller BIC
            
            min_i <- model_i[ which.min( BIC_vec[model_i] ) ]
            
            BIC_vec_new[i] <- BIC_vec[min_i]
            AIC_vec_new[i] <- AIC_vec[min_i]
            fit_list_new[[i]] <- fit_list[[min_i]]
            mu_list_new[[i]] <- mu_list[[min_i]]
            pi_list_new[[i]] <- pi_list[[min_i]]
            pi0_list_new[[i]] <- pi0_list[[min_i]]
            gamma_list_new[[i]] <- gamma_list[[min_i]]
        }
    }    
    
    return( list( fit_list=fit_list_new, mu_list=mu_list_new,
        pi_list=pi_list_new, pi0_list=pi0_list_new, 
        gamma_list=gamma_list_new,
        BIC_vec=BIC_vec_new, AIC_vec=AIC_vec_new ) )
}
