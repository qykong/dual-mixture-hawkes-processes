# This script implements the dual mixture variant of the power law decayed Hawkes processes: a kernel mixture model and a Borel mixture model
library(tidyverse)
library(evently)
library(parallel)
source('scripts/utility-functions.R')
source('scripts/kmm_m_step.R')

######### kernel mixture model fitting ##################
KMMEM <- function(data, k, max_iter = 10) {
  ps <- random_init_probabilities(k)
  
  params <- generate_random_points.hawkes_mixturePLkernel(new_hawkes(model_type='mixturePLkernel', model_vars = list(cluster_number = k)))[1,]
  params <- params %>% as.double() %>% set_names(names(params))
  
  data <- keep(data, ~nrow(.x) > 1)
  kernel_log_likelihood <- function(cluster_no, param, hist) {
    c <- param[[sprintf('c[%s]', cluster_no)]]
    theta <- param[[sprintf('theta[%s]', cluster_no)]]
    if (nrow(hist) == 1) stop()
    sum(map_dbl(seq(2, nrow(hist)), ~log(sum(theta * c^theta * (hist$time[.x] - hist$time[seq(.x-1)] + c) ^ (-1-theta)) +.Machine$double.xmin )))
  }

  iter <- 1
  total_likelihood <- -.Machine$double.xmax
  print('start em on kernel functions....')
  
  repeat {
    qs <- map(seq_along(data), function(l) map_dbl(seq(k), function(.x) (log(ps[.x]) + kernel_log_likelihood(.x, params, data[[l]]))))
    new_total_likelihood <- sum(map_dbl(seq_along(data), function(l) log_sum_exp(qs[[l]]) ))
    print(sprintf('iteration %s: %s', iter, round(new_total_likelihood, digits = 2)))
    if (abs(new_total_likelihood - total_likelihood) < 1e-2 || iter >= max_iter) {
      break
    } else {
      iter <- iter + 1
      total_likelihood <- new_total_likelihood
    }
    p_k_l <- map(seq(k), function(.x) map_dbl(seq_along(data), function(l) if (is.infinite(qs[[l]][.x])) 0 else 1/(sum(exp(qs[[l]][-.x] - qs[[l]][.x])) + 1) ))
    p_k_h_index <- map_df(seq(k), function(.x) {
      map_dfr(seq_along(data), function(l) list(k = .x, h = l, p = p_k_l[[.x]][[l]]))
    })
    if (any(is.nan(p_k_h_index$p))) browser()
    ps <- map_dbl(seq(k), function(.x) sum(p_k_l[[.x]])/length(data))
    res <- fit_series(data, model_type = 'mixturePLkernel', init_pars = set_names(rbind.data.frame(params), names(params)),
                      observation_time = Inf, model_vars = list(cluster_number = k, p_k_h_index = p_k_h_index), cores = 1)
    params <- res$par
  }
  params <- map(seq(k), ~c(theta = params[[sprintf('theta[%s]', .x)]],
                           c = params[[sprintf('c[%s]', .x)]]))
  return(list(probability = ps, params = params, p_k_h_index = p_k_h_index, total_likelihood = total_likelihood))
}

KMMEM_repeat <- function(..., times = 10, cores = 1) {
  models <- mclapply(seq(times), function(i) {
    tryCatch({
      KMMEM(...)
    }, error = function(e) {
      print(e)
      return(list(total_likelihood = -.Machine$double.xmax))
    })
  }, mc.cores = cores)
  models[[which.max(map_dbl(models, 'total_likelihood'))]]
}


######### Borel mixture model fitting ##################

BMMEM <- function(sizes, k = 1, max_iter = 50) {
  counts <- as.data.frame(table(sizes), stringsAsFactors = F) %>%
    rename(size = sizes, count = Freq) %>%
    mutate(size = as.integer(size))
  ps <- random_init_probabilities(k)
  n_stars <- runif(k, 0, 1)
  
  borel <- function(x, lam) {
    dpois(x-1, x*lam)/x
  }
  iter <- 1
  total_likelihood <- -.Machine$double.xmax
  repeat {
    qs <- map(seq(k), function(.x) map_dbl(seq_along(counts$size), function(l) ps[.x] * borel(counts$size[l], n_stars[.x])))
    new_total_likelihood <- sum(map_dbl(seq_along(counts$size), function(l) log(sum(map_dbl(seq(k), function(.x) qs[[.x]][l])))) * counts$count)
    if (new_total_likelihood - total_likelihood < 1e-2 || iter >= max_iter) {
      break
    } else {
      iter <- iter + 1
      total_likelihood <- new_total_likelihood
    }
    
    s_q <- map_dbl(seq_along(counts$size), function(l) sum(map_dbl(seq(k), function(.x) qs[[.x]][l])))
    p_k_l <- map(seq(k), function(.x) map_dbl(seq_along(counts$size), function(l) qs[[.x]][l]/s_q[l] ))
    n_stars <- map_dbl(seq(k), function(.x) sum((counts$size-1)*p_k_l[[.x]] * counts$count)/sum((counts$size)*p_k_l[[.x]] * counts$count) )

    ps <- map_dbl(seq(k), function(.x) sum(p_k_l[[.x]] * counts$count)/sum(counts$count) )
  }
  
  return(list(n_star = n_stars,
              p = ps,
              total_likelihood = total_likelihood))
}

BMMEM_repeat <- function(..., times = 10, cores = 1) {
  models <- mclapply(seq(times), function(i) {BMMEM(...) }, mc.cores = cores)
  models[[which.max(map_dbl(models, 'total_likelihood'))]]
}

######## fit the two mixture models together #########
clustering_n_star_and_kernel <- function(hists, clusters, cores = 1, times = 10) {
  sizes <- map_int(hists, nrow)
  cat('start BMM\n')
  if (!missing(clusters)) {
    BMM_clusters <- clusters[1]
    n_star_p <- BMMEM_repeat(sizes, k = BMM_clusters, cores = cores, times = times)
    if (length(clusters) == 2) kernel_clusters <- clusters[2] else kernel_clusters <- clusters
  } else {
    # test k from 2 to 10 and choose the best by AIC
    trials <- map(seq(2, 10), function(k) {
      n_star_p <- BMMEM_repeat(sizes, k = k, cores = cores, times = times)
      list(n_star_p = n_star_p,
           aic_n_star = -n_star_p$total_likelihood * 2 + k * 2,
           k = k)
    })
    best_ind <- which.min(map_dbl(trials, 'aic_n_star'))
    BMM_clusters <- trials[[best_ind]]$k
    kernel_clusters <- trials[[best_ind]]$k
    n_star_p <- trials[[best_ind]]$n_star_p
    cat(sprintf('best number of clusters is %s\n', BMM_clusters))
  }
  
  cat('BMM done; start clustering kernel functions\n')
  
  # remove single event cascades as they won't be computed in KMM anyway
  keeped_hists <- keep(hists, ~nrow(.x) >= 2)
  # if no cascades left then return here
  if (length(keeped_hists) == 0) return(list(n_star = n_star_p$n_star,
                                             p_n_star = n_star_p$p,
                                             params = NA,
                                             p_params = NA,
                                             kernel_clusters = kernel_clusters,
                                             BMM_clusters = BMM_clusters))
  cat('start KMM\n')
  kernel_clusters <- min(length(keeped_hists), kernel_clusters)
  
  res <- KMMEM_repeat(keeped_hists, k = kernel_clusters, times = times, cores = cores)
  cat('done KMM\n')
  params <- res
  return(list(n_star = n_star_p$n_star,
              p_n_star = n_star_p$p,
              params = res$params,
              p_params = res$probability,
              kernel_clusters = kernel_clusters,
              BMM_clusters = BMM_clusters))
}