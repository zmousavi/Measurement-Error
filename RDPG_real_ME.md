```
---
title: "RDPG_ME_size"
output:
  md_document:
    variant: markdown_github
---
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


Real RDPG - 100 mc runs
red: naive regression coefficient estimate
green: adjusted regression coefficient estimate
blue: true regression coefficient estimate

beta0 coefficient:


![](RDPG_ME_size_files/figure-markdown_github/beta0_size.jpg)

beta1 coefficient:
![](RDPG_ME_size_files/figure-markdown_github/beta1_size.jpg)

MSE Ratio (all graph sizes):
![](RDPG_ME_size_files/figure-markdown_github/beta_real_strc_mse_size.jpg)

MSE Ratio (graphs larger than 500:
![](RDPG_ME_size_files/figure-markdown_github/beta_real_strc_mse_size_500.jpg)

delta_bias:
![](RDPG_ME_size_files/figure-markdown_github/delta_bias_real_strc_size.jpg)

delta_variance:
![](RDPG_ME_size_files/figure-markdown_github/delta_var_real_strc_size.jpg)

mean(|Xhat1 - X1|):
![](RDPG_ME_size_files/figure-markdown_github/Xhat1_bias_size.jpg)


