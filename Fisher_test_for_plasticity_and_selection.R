#low latitude populations

#proportions
dat1 <- matrix(c(0.338,0.204))
#o fisher nao da para fazer
fisher.test(dat1)
# o chi-square da
chisq.test(dat1)

#other way but diferent from what we did in the paper
dat2 <- data.frame(
  "Selection_yes" = c(321, 628),
  "Selection_no" = c(2106, 8845),
  row.names = c("Plastic_yes", "Plastic_no"),
  stringsAsFactors = FALSE
)
fisher.test(dat2)
