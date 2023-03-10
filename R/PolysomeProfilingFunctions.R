normalize <- function(dtf, max_abs = Inf, pos_start = 7,
                      pos_end = 60, pos_offset = 0,
                      to = 'AUC', max_jump = Inf,
                      smoothen = TRUE){

  dtf2 <- dtf %>%
    filter(`Position(mm)` > pos_start) %>%
    filter(`Position(mm)` < pos_end) %>%
    mutate(`Position(mm)` = `Position(mm)` + pos_offset) %>%
    mutate(Absorbance = Absorbance + abs(min(Absorbance)))

  dtAbs <- abs(dtf2$Absorbance[1:length(dtf2$Absorbance)-1] - dtf2$Absorbance[2:length(dtf2$Absorbance)])
  dtf2$dAbs <- append(dtAbs, 0)
  dtf2 <- dtf2 %>%
    filter(dtf2$dAbs < max_jump)

  if(smoothen){
    fit.abs <- smooth.spline(x = dtf2$`Position(mm)`, y = dtf2$Absorbance, df = 100)
    dtf2$Absorbance <- predict(fit.abs, x = dtf2$`Position(mm)`)$y
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
                  ref_start = 270, ref_end = 940, by_peaks = TRUE, npeaks = 4, show_peaks = FALSE, minPeakPos = 15, maxPeakPos = 70){
  if(by_peaks){
    other$`Position(mm)` <- other$`Position(mm)` - mean(peak_finder(other, npeaks, show_peaks, minPeakPos, maxPeakPos)$`Position(mm)` - peak_finder(ref, npeaks, show_peaks, minPeakPos, maxPeakPos)$`Position(mm)`)
  }else{
    other$`Position(mm)` <- other$`Position(mm)` - calc_offset(ref, other, span = span, return_plot=return_plot, ref_start = ref_start, ref_end = ref_end) * 0.061
  }
  return(other)
}

PolyMonoRatio <- function(dtf, mono_start=29, mono_end=37, poly_start=mono_end, poly_end=60, show_cutoffs=TRUE){
  dtf_mono <- dtf %>%
    filter(`Position(mm)` > mono_start & `Position(mm)` < mono_end)
  dtf_poly <- dtf %>%
    filter(`Position(mm)` > poly_start & `Position(mm)` < poly_end)

  if(show_cutoffs){
    p <- ggplot(dtf, mapping = aes(x = `Position(mm)`, y = Absorbance)) +
      geom_line(size = 0.6) +
      geom_hline(yintercept = 0, color = 'black') +
      geom_vline(xintercept = mono_start, color = 'red') +
      geom_vline(xintercept = mono_end, color = 'red') +
      geom_vline(xintercept = poly_start, color = 'blue') +
      geom_vline(xintercept = poly_end, color = 'blue') +
      theme_bw()

    print(p)
  }

  return(trapz(dtf_poly$`Position(mm)`, dtf_poly$Absorbance) /
           trapz(dtf_mono$`Position(mm)`, dtf_mono$Absorbance))
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


PrismExport <- function(dtf, wider_names = c('Sample_ID'), wider_vals = c('Absorbance'), mode = 'percentage'){

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

  return(to_exp)
}


peak_finder <- function(dtf, npeaks = 4, show_peaks = FALSE, minPeakPos = 15, maxPeakPos = 70){
  tst1 <- dtf
  fit.abs <- smooth.spline(x = tst1$`Position(mm)`, y = tst1$Absorbance, df = 50)
  tst1$Absorbance <- predict(fit.abs, x = tst1$`Position(mm)`)$y

  tst1$Absorbance <- tst1$Absorbance + abs(min(tst1$Absorbance))
  peaks <- findpeaks(tst1$Absorbance, peakpat = '[+]{30,}[-]{30,}')
  peaks <- data.frame(peaks)
  peaks$pos <- tst1$`Position(mm)`[peaks$X2]

  good_peaks <- peaks %>%
    filter(pos > minPeakPos, pos < maxPeakPos) %>%
    rename(absorb = X1) %>%
    slice_head(n = npeaks)

  peaks$good <- do.call(paste0, peaks) %in% do.call(paste0, good_peaks)

  if(show_peaks){
    p <- ggplot(data = tst1, aes(x = `Position(mm)`, y = Absorbance)) +
      geom_line() +
      geom_point(data = peaks, aes(x = pos, y = X1, color = good), size = 3) +
      theme_bw()
    print(p)
  }
  return(dtf[good_peaks$X2, ])
}

