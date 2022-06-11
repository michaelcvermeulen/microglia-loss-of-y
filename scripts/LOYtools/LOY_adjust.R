LOY_adjust <- function(dat, outliers=100, plot = F, min_cells = 100){
  
  library(aomisc)
  readr::read_csv(
  "https://raw.githubusercontent.com/michaelcvermeulen/microglia-loss-of-y/main/data/LOY_tables/all_cell_types/LOY_prop_table_3000_1000_.txt") %>%
   as.data.table() -> l

  filter_function(min_cells = min_cells, score = -10, l = l, type = l$cell_type %>% unique()) -> d
  
  
    model <- drc::drm(LOY_percent ~ sum_Y_exp, fct = aomisc::DRC.expoDecay(),
             data = d[d$LOY_percent<=outliers,])
    predict(model, newdata = dat) -> dat$predicted
    
    dplyr::mutate(dat, residuals = LOY_percent - predicted) -> dat
    ifelse(dat$residuals<0, 0 , dat$residuals) -> dat$adj_LOY_percent
   
    

  if(plot==T){ 
  
    plot(model, log = "")
  
    dat %>%
    dplyr::mutate(curve = predict(model, newdata = dat)) %>%
    ggplot(aes(sum_Y_exp,LOY_percent)) +
    geom_point(color = "grey50") +
    geom_line(aes(y = curve)) -> a 
    
    a
  }

    return(list(dat,model))
}


# function to pull p value from lm model
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}