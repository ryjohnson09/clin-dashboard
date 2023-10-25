efficacy_models <- function(data, var=NULL, wk=NULL, show_pvalue = TRUE) {
  # Need to set contrasts to work for Type III SS. See analysis results metadata for
  # table 14-3.01. Reference for R here: https://www.r-bloggers.com/anova-%E2%80%93-type-iiiiii-ss-explained/
  op <- options(contrasts = c("contr.sum","contr.poly"))
  
  # Subset to analyze
  data <- data %>%
    filter(AVISITN == wk)
  
  data <- data %>%
    mutate(
      TRTPCD = case_when(
        TRTPN == 0 ~ 'Pbo',
        TRTPN == 54 ~ 'Xan_Lo',
        TRTPN == 81 ~ 'Xan_Hi'
      )
    )
  
  # Create an ordered factor variable for the models
  data['TRTPCD_F'] <- factor(data$TRTPCD, levels=c('Xan_Hi', 'Xan_Lo', 'Pbo'))
  data['AWEEKC'] = factor(data$AVISIT)
  
  # Set up the models
  if (var == "CHG") {
    model1 <- lm(CHG ~ TRTPN + SITEGR1 + BASE, data=data)
    model2 <- lm(CHG ~ TRTPCD_F + SITEGR1 + BASE, data=data)
  } else {
    model1 <- lm(AVAL ~ TRTPN + SITEGR1, data=data)
    model2 <- lm(AVAL ~ TRTPCD_F + SITEGR1, data=data)
  }
  
  ## Dose Response --- NOTE: For statistics portions, I purposefully did not
  #import the libraries to make it explicitly clear which packages were being
  #used to match P-values.
  ancova <- drop1(model1, .~., test="F")
  
  # Pull it out into a table
  sect1 <- tibble::tibble(row_label=c('p-value(Dose Response) [1][2]'),
                          `81` = ifelse(show_pvalue, c(num_fmt(ancova[2, 'Pr(>F)'], int_len=4, digits=3, size=12)), "Not Applicable")
  ) %>%
    pad_row()
  
  ## Pairwise Comparisons ----
  # Here's a reference for the emmeans package and how to use it:
  #   https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html
  # Adjustments made are in line with the analysis results metadata in the analysis define
  # and PROC GLM documentation.
  
  # Linear model but use treatment group as a factor now
  # LS Means and weight proportionately to match OM option on PROC GLM in SAS
  lsm <- emmeans::lsmeans(model2, ~TRTPCD_F, weights='proportional')
  
  # Here on out - it's all the same data manipulation
  # Get pairwise contrast and remove P-values adjustment for multiple groups
  cntrst_p <- emmeans::contrast(lsm, method="pairwise", adjust=NULL)
  # 95% CI
  cntrst_ci <- confint(cntrst_p)
  
  # merge and convert into dataframe
  pw_data <- tibble::as_tibble(summary(cntrst_p)) %>%
    merge(tibble::as_tibble(cntrst_ci)) %>%
    rowwise() %>%
    # Create the display strings
    mutate(
      p = ifelse(show_pvalue, num_fmt(p.value, int_len=4, digits=3, size=12), "Not Applicable"),
      diff_se = as.character(
        glue::glue('{num_fmt(estimate, int_len=2, digits=1, size=4)} ({num_fmt(SE, int_len=1, digits=2, size=4)})')
      ),
      ci = as.character(
        glue::glue('({num_fmt(lower.CL, int_len=2, digits=1, size=4)};{num_fmt(upper.CL, int_len=1, digits=1, size=3)})')
      )
    ) %>%
    # Clean out the numeric variables
    select(contrast, p, diff_se, ci) %>%
    # Transpose
    tidyr::pivot_longer(c('p', 'diff_se', 'ci'), names_to = 'row_label')
  
  # Subset Xan_Lo - Pbo into table variables
  xan_lo <- pw_data %>%
    filter(contrast == 'Xan_Lo - Pbo') %>%
    # Rename to the table display variable
    select(`54`=value) %>%
    pad_row()
  
  #Add in row_label
  xan_lo['row_label'] <- c('p-value(Xan - Placebo) [1][3]', '  Diff of LS Means (SE)', '  95% CI', '')
  
  # Subset Xan_hi - Pbo into table variables
  xan_hi <- pw_data %>%
    filter(contrast == 'Xan_Hi - Pbo') %>%
    # Rename to the table display variable
    select(`81`=value) %>%
    pad_row()
  # Add in row_label
  xan_hi['row_label'] <- c('p-value(Xan - Placebo) [1][3]', '  Diff of LS Means (SE)', '  95% CI', '')
  xan_hi['ord'] <- c(1,2,3,4) # Order for sorting
  
  # Subset Xan_Hi - Xan_Lo into table variable
  xan_xan <- pw_data %>%
    filter(contrast == 'Xan_Hi - Xan_Lo') %>%
    # Rename to the table display variable
    select(`81`=value)
  # Add in row_label
  xan_xan['row_label'] <- c('p-value(Xan High - Xan Low) [1][3]', '  Diff of LS Means (SE)', '  95% CI')
  xan_xan['ord'] <- c(5,6,7) # Order for sorting
  
  # Pack it all together
  pw_final <- merge(xan_lo, xan_hi, by='row_label') %>%
    bind_rows(xan_xan) %>%
    arrange(ord)
  
  # Bind and clean up
  bind_rows(sect1, pw_final) %>% 
    select(row_label, 
           `var1_Xanomeline Low Dose` = `54`,
           `var1_Xanomeline High Dose` = `81`
    )
}

#####################################

pad_row <- function(.data, n=1) {
  .data[(nrow(.data)+1):(nrow(.data)+n), ] <- ""
  .data
}

#####################################

nest_rowlabels <- function(.dat) {
  stubs <- .dat %>% 
    distinct(row_label1, ord_layer_index) %>% 
    rename(row_label = row_label1) %>% 
    mutate(
      ord_layer_1 = 0, 
      ord_layer_2 = 0
    )
  
  .dat %>% 
    select(-row_label1, row_label=row_label2) %>% 
    bind_rows(stubs) %>% 
    arrange(ord_layer_index, ord_layer_1, ord_layer_2) %>% 
    mutate(
      across(starts_with('var'), ~ tidyr::replace_na(., ''))
    )
}

#####################################

num_fmt <- Vectorize(function(var, digits=0, size=10, int_len=3) {
  # Formats summary stat strings to align display correctly
  
  if (is.na(var)) return('')
  
  # Set nsmall to input digits
  nsmall = digits
  
  # Incremement digits for to compensate for display
  if (digits > 0) {
    digits = digits + 1
  }
  
  # Form the string
  return(str_pad(
    format(
      # Round
      round(var, nsmall),
      # Set width of format string
      width=(int_len+digits),
      # Decimals to display
      nsmall=nsmall
    ),
    # Overall width padding
    side='right', size
  ))
})

##########################################

fmt_est <- function(.mean,
                    .sd,
                    digits = c(1, 2)) {
  .mean <- fmt_num(.mean, digits[1], width = digits[1] + 4)
  .sd <- fmt_num(.sd, digits[2], width = digits[2] + 3)
  paste0(.mean, " (", .sd, ")")
}

##########################################

fmt_ci <- function(.est,
                   .lower,
                   .upper,
                   digits = 2,
                   width = digits + 3) {
  .est <- fmt_num(.est, digits, width)
  .lower <- fmt_num(.lower, digits, width)
  .upper <- fmt_num(.upper, digits, width)
  paste0(.est, " (", .lower, ",", .upper, ")")
}

##########################################

fmt_num <- function(x, digits, width = digits + 4) {
  formatC(x,
          digits = digits,
          format = "f",
          width = width
  )
}

##########################################

fmt_pval <- function(.p, digits = 3) {
  scale <- 10^(-1 * digits)
  p_scale <- paste0("<", digits)
  ifelse(.p < scale, p_scale, fmt_num(.p, digits = digits))
}