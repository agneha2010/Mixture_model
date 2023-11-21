options(kableExtra.latex.load_packages = FALSE) 
library(kableExtra)
kbl(dt,"latex")
kbl(dt, booktabs = T, linesep = "","latex") %>% 
  kable_styling(latex_options = c("striped")) %>% 
  add_header_above(c("semi parametric" = 1, "parametric" = 3, "semi parametric" = 1, "parametric" = 3)) %>% 
  add_header_above(c("beta" = 4, "a" = 4)) %>% 
  kable_styling(position = "center") %>% 
  column_spec(5, border_left = T)

kbl(dt, booktabs = T, linesep = "") %>% 
  kable_styling(latex_options = c("striped")) %>% 
  add_header_above(c("semi parametric" = 1, "parametric" = 3, "semi parametric" = 1, "parametric" = 3)) %>% 
  add_header_above(c("beta" = 4, "a" = 4)) %>% 
  kable_styling(position = "center") %>% 
  column_spec(5, border_left = T)
