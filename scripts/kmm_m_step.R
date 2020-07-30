# this script implements the M step of learning KMM (power law decay) with evently (optimized by AMPL)

# hack due to a bug in evently (get_param_names and get_ampl_likelihood are not exposed)
get_param_names <- getFromNamespace('get_param_names', ns = 'evently')
get_ampl_likelihood <- getFromNamespace('get_ampl_likelihood', ns = 'evently')

# defines the log-likelihood function
get_ampl_likelihood.hawkes_mixturePLkernel <- function(model) {
  paste(
    'sum {cn in 1..HL} ( sum {cli in 1..CL} (  PKH[cli, cn] * (1 - 0^(L[cn]-J0[cn]-1)) *  (sum {i in J0[cn]+1..L[cn]-1} (',
    'log(sum {j in 1..i-1} (theta[cli] * c[cli]^theta[cli] * (time[cn,i] - time[cn,j] + c[cli]) ^ (-1-theta[cli])) + 1e-100)',
    '  ))',
    '));'
  )
}

# customized AMPL data output
get_ampl_data_output.hawkes_mixturePLkernel <- function(model) {
  output <- NextMethod()
  c(output,
    evently:::ampl_output_from_r('PKH', model$p_k_h_index, 'data.frame'))
}

# customized AMPL model output
get_ampl_model_output.hawkes_mixturePLkernel <- function(model) {
  paste('param HL > 0; param CL := ', model$cluster_number, '; param ML > 0; param L {1..HL} >= 0; param magnitude {1..HL,1..ML} >= 0; param time {1..HL,1..ML} >= 0; param J0 {1..HL} >= 0; param PKH {1..CL,1..HL} >=0;
        var c {1..CL} >= 0;
        var theta {1..CL} >= 0;',
        if (is.finite(model$init_par[[1]])) paste(glue::glue('let {get_param_names(model)} := {model$init_par[get_param_names(model)]};'), collapse = "\n") else '',
        'maximize Likelihood:',
        get_ampl_likelihood(model),
        # parameter constraints
        paste(glue('subject to c_limit_{seq(model$cluster_number)}: c[{seq(model$cluster_number)}] <= 200;'), collapse = '\n'),
        paste(glue('subject to theta_limit_{seq(model$cluster_number)}: theta[{seq(model$cluster_number)}] <= 5;'), collapse = '\n'),
        sep = '\n'
  )
}

# defines the parameters (`model$cluster_number` sets of theta and c)
get_param_names.hawkes_mixturePLkernel <- function(model) {
  c(sprintf('theta[%s]', seq(model$cluster_number)),
    sprintf('c[%s]', seq(model$cluster_number)))
}

# provide random initialisations for c and theta
generate_random_points.hawkes_mixturePLkernel <- function(model) {
  params <- get_param_names(model)
  inits <- map(seq_along(params), ~runif(10, min = .Machine$double.eps, max = 3)) %>%
    {do.call(cbind, .)} %>%
    as.data.frame() %>%
    set_names(params)
  inits[9, ] <- NA
  inits[10, ] <- Inf
  inits
}
