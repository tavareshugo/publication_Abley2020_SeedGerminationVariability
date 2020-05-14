# Function to sample null distributions of days
#' @param germ_day a vector of germination days
#' @param prob_day a vector of probabilities of each germination day
#' @param bottom_n,top_n integer indicating the number of seeds from each half of the silique
sample_germ <- function(germ_day, prob_day, bottom_n, top_n){
  
  # Check input parameters
  if(length(bottom_n) != 1 | length(top_n) != 1) stop("bottom_n and top_n should be of length 1.")
  if(bottom_n %% 1 != 0 | top_n %% 1 != 0) stop("bottom_n and top_n need to be integers.")
  if(length(germ_day) != length(prob_day)) stop("germ_day and prob_day have to be of same length.")
  if(any(prob_day < 0 | prob_day > 1)) stop("prob_day should only contain values between 0 and 1.")
  if(any(germ_day < 1) | any(germ_day %% 1 != 0)) stop("germ_day should only contain integers greater than zero.")
  
  # Sample germination days for bottom silique
  bottom_sample <- sample(germ_day, bottom_n, prob = prob_day, replace = TRUE)
  
  # Sample germination days for top silique
  top_sample <- sample(germ_day, top_n, prob = prob_day, replace = TRUE)
  
  # output a table
  data.frame(sim_cv_bottom = sd(bottom_sample)/mean(bottom_sample),
             sim_cv_top = sd(top_sample)/mean(top_sample),
             sim_d20_bottom = sum(bottom_sample >= 20)/bottom_n,
             sim_d20_top = sum(top_sample >= 20)/top_n)
}

# Function to bootstrap this sampling
boot_germ <- function(germ_day, prob_day, bottom_n, top_n, nboot = 1000){
  purrr::rerun(nboot,
               sample_germ(germ_day, prob_day, bottom_n, top_n)) %>% 
    bind_rows()
}

# # Test the function
# boot_germ(1:5, c(0.1, 0.2, 0.1, 0.2, 0.4), 5, 20) %>%
#   ggplot() + geom_histogram(aes(cv_bottom, fill = "bottom"), alpha = 0.5) +
#   geom_histogram(aes(cv_top, fill = "top"), alpha = 0.5)
