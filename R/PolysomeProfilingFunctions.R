normalize <- function(dtf, max_abs = Inf, pos_start = 7,
                      pos_end = 67, pos_offset = 0,
                      to = 'AUC', max_jump = Inf,
                      smoothen = TRUE, zero_baseline = TRUE){

  dtf2 <- dtf %>%
    filter(`Position(mm)` > pos_start) %>%
    filter(`Position(mm)` < pos_end) %>%
    mutate(`Position(mm)` = `Position(mm)` + pos_offset)

  dtAbs <- abs(dtf2$Absorbance[1:length(dtf2$Absorbance)-1] - dtf2$Absorbance[2:length(dtf2$Absorbance)])
  dtf2$dAbs <- append(dtAbs, 0)
  dtf2 <- dtf2 %>%
    filter(dtf2$dAbs < max_jump) %>%
    mutate(Absorbance = Absorbance - min(Absorbance) * zero_baseline)

  if(smoothen){
    fit.abs <- smooth.spline(x = dtf2$`Position(mm)`, y = dtf2$Absorbance, df = 100)
    dtf2$Absorbance <- predict(fit.abs, x = dtf2$`Position(mm)`)$y
    dtf2 <- dtf2 %>% mutate(Absorbance = Absorbance - min(Absorbance) * zero_baseline)
  }

  if(to == '80S'){
    temp <- dtf2 %>%
      dplyr::filter(`Position(mm)` > 20)

    dtf2 <- dtf2 %>%
      mutate(Absorbance = Absorbance/max(temp$Absorbance))
  }

  if(to == 'AUC'){
    dtf2 <- dtf2 %>%
      mutate(Absorbance = Absorbance/trapz(dtf2$`Position(mm)`, dtf2$Absorbance))
  }

  if(to == 'AUC_Rib'){
    temp <- dtf2 %>%
      dplyr::filter(`Position(mm)` > 20)

    dtf2 <- dtf2 %>%
      mutate(Absorbance = Absorbance/trapz(temp$`Position(mm)`, temp$Absorbance))
  }

  dtf2 <- dtf2 %>%
    mutate(Absorbance = replace(Absorbance, Absorbance > max_abs, max_abs))
  return(dtf2)
}

calc_offset <- function(ref, other, span = seq(-250, 250, 1), return_plot=FALSE,
                        ref_start = 270, ref_end = 940){
  dtfAbsDiff <- c()
  for(i in span){
    dtfAbsDiff <- append(dtfAbsDiff, sum((other$Absorbance[(ref_start+i):(ref_end+i)] - ref$Absorbance[ref_start:ref_end])^2))
  }
  if(return_plot){
    plot(span, dtfAbsDiff)
  }
  return(span[which.min(dtfAbsDiff)])
}

align <- function(ref, other, span = seq(-250, 250, 1), return_plot=FALSE,
                  ref_start = 270, ref_end = 940, by_peaks = TRUE, npeaks = 1, show_peaks = FALSE, minPeakPos = 26, maxPeakPos = 46){
  if(by_peaks){
    other$`Position(mm)` <- other$`Position(mm)` - mean(peak_finder(other, npeaks, show_peaks, minPeakPos, maxPeakPos)$`Position(mm)` - peak_finder(ref, npeaks, show_peaks, minPeakPos, maxPeakPos)$`Position(mm)`)
  }else{
    other$`Position(mm)` <- other$`Position(mm)` - calc_offset(ref, other, span = span, return_plot=return_plot, ref_start = ref_start, ref_end = ref_end) * 0.061
  }
  return(other)
}

PolyMonoRatio <- function(dtf, mono_start=29, mono_end=37, poly_start=mono_end, poly_end=60, show_cutoffs=TRUE, adjust_baseline=FALSE){

  if(show_cutoffs){
    p <- ggplot(dtf, mapping = aes(x = `Position(mm)`, y = Absorbance)) +
      geom_line(size = 0.6) +
      geom_vline(xintercept = mono_start, color = 'red') +
      geom_vline(xintercept = mono_end, color = 'red') +
      geom_vline(xintercept = poly_start, color = 'blue') +
      geom_vline(xintercept = poly_end, color = 'blue') +
      theme_bw()
    }

  if(adjust_baseline){

    monos <- peak_area(dtf, start_pos = mono_start, end_pos = mono_end)
    polys <- peak_area(dtf, start_pos = poly_start, end_pos = poly_end)

    if(show_cutoffs){
      print(p + monos[2] + polys[2])
    }

    return(polys[[1]] / monos[[1]])
  }

  dtf_mono <- dtf %>%
    filter(`Position(mm)` > mono_start & `Position(mm)` < mono_end)
  dtf_poly <- dtf %>%
    filter(`Position(mm)` > poly_start & `Position(mm)` < poly_end)

  if(show_cutoffs){
    print(p + geom_hline(yintercept = 0, color = 'black'))
  }

  return(trapz(dtf_poly$`Position(mm)`, dtf_poly$Absorbance) /
           trapz(dtf_mono$`Position(mm)`, dtf_mono$Absorbance))
}


peak_area <- function(dtf, start_pos=29, end_pos=37){

  poly_abs <- subset(dtf, `Position(mm)` > start_pos & `Position(mm)` < end_pos)
  fit_poly <- subset(poly_abs, Absorbance == min(Absorbance[1:length(Absorbance)%/%5]) |
                       Absorbance == min(Absorbance[(length(Absorbance)%/%2 + 1):length(Absorbance)]))
  fit_poly_abs <- lm(Absorbance ~`Position(mm)`, data = fit_poly)

  poly_abs$Absorbance2 <- poly_abs$Absorbance - predict(fit_poly_abs, newdata = poly_abs)
  poly_abs$Absorbance2[poly_abs$Absorbance2 < 0] <- 0

  prof_baseline <- geom_line(data = poly_abs, aes(x = `Position(mm)`, y = predict(fit_poly_abs, newdata = poly_abs)), color = 'green')

  return(c(trapz(poly_abs$`Position(mm)`, poly_abs$Absorbance2), prof_baseline))
}


ProfLoader <- function(excel_path, excel_sheet, excel_range, prof_path, prof_pattern, output_as = 'list'){

  dtf1 <- read_xlsx(excel_path, sheet = excel_sheet, range = excel_range)

  profList <- list.files(path = prof_path, pattern = prof_pattern)
  dtfs_prof <- lapply(profList, function(x) {read_csv(paste(prof_path, x, sep = '/'), skip = 47)})
  for(i in c(1:length(dtfs_prof))){
    dtfs_prof[[i]]$file_name <- profList[i]
  }

  all_raw1 <- left_join(dtf1, rbindlist(dtfs_prof), by = c('Sample_name' = 'file_name'))

  if(output_as == 'list'){
    return(group_split(all_raw1, Sample_ID))
  }
  return(all_raw1)
}


PrismExport2 <- function(dtf, wider_names = c('Sample_ID'), wider_vals = c('Absorbance'), mode = 'percentage'){

  to_exp <- dtf

  if(mode == 'percentage'){
    to_exp <- lapply(to_exp, function(x) mutate(x,
                                                Absorbance = Absorbance + abs(min(Absorbance)),
                                                Absorbance = Absorbance/max(Absorbance)*100))
  }

  to_exp <- rbindlist(to_exp) %>%
    select(c('Position(mm)', wider_names, wider_vals)) %>%
    pivot_wider(names_from = wider_names, values_from = wider_vals) %>%
    arrange(`Position(mm)`)

  for (v in colnames(to_exp)){
    to_exp[v] <- as.numeric(as.character(unlist(to_exp[v], recursive = FALSE)))
  }

  return(to_exp)
}


peak_finder <- function(dtf, npeaks = 1, show_peaks = FALSE, minPeakPos = 26, maxPeakPos = 46, minAbs = 0){
  tst1 <- dtf
  fit.abs <- smooth.spline(x = tst1$`Position(mm)`, y = tst1$Absorbance, df = 100)
  tst1$Absorbance <- predict(fit.abs, x = tst1$`Position(mm)`)$y

  tst1$Absorbance <- tst1$Absorbance + abs(min(tst1$Absorbance))
  peaks <- findpeaks(tst1$Absorbance, peakpat = '[+]{30,}[-]{30,}')
  peaks <- data.frame(peaks)
  peaks$pos <- tst1$`Position(mm)`[peaks$X2]

  good_peaks <- peaks %>%
    filter(pos > minPeakPos, pos < maxPeakPos, X1 > minAbs) %>%
    rename(absorb = X1) %>%
    slice_head(n = npeaks)

  peaks$good <- do.call(paste0, peaks) %in% do.call(paste0, good_peaks)

  if(show_peaks){
    p <- ggplot(data = tst1, aes(x = `Position(mm)`, y = Absorbance)) +
      geom_line() +
      geom_point(data = peaks, aes(x = pos, y = X1, color = good), size = 3) +
      geom_text(data = peaks, aes(x = pos, y = X1),
                aes(label = paste(round(pos,2), round(X1, 2), sep = ', '))) +
      theme_bw()
    print(p)
  }
  return(dtf[good_peaks$X2, ])
}

