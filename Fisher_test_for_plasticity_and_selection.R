#no input file is needed
#output: gives

#low latitude populations

dat1 <- data.frame(
  "Selection_yes" = c(321, 628),
  "Selection_no" = c(2106, 8845),
  row.names = c("Plastic_yes", "Plastic_no"),
  stringsAsFactors = FALSE
)
fisher.test(dat1)


#high latitude populations

dat2 <- data.frame(
  "Selection_yes" = c(723, 3472),
  "Selection_no" = c(844, 6861),
  row.names = c("Plastic_yes", "Plastic_no"),
  stringsAsFactors = FALSE
)
fisher.test(dat2)
