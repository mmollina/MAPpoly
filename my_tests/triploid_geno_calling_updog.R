require(vcfR)
require(updog)
require(tidyverse)
require(VariantAnnotation)
dat <- vcfR::read.vcfR("~/repos/official_repos/MAPpoly/my_tests/BreedBaseGenotypesDownload.vcf")
ad <- extract.gt(dat, element = "AD")
ad.f1 <- ad[,str_detect(colnames(ad), pattern = "F1")]
ad.f1.1 <- masplit(ad, record = 1)
ad.f1.2 <- masplit(ad, record = 2)
res <- multidog(refmat = ad.f1.1, 
                sizemat = ad.f1.1 + ad.f1.2, 
                model = "hw", 
                ploidy = 3, 
                nc = 32)
plot(res, indices = 2)
saveRDS(res, file = "~/repos/official_repos/MAPpoly/my_tests/updog_banana.rds")
## Monyet: 4X
## Kokopo: 2X
ad[2,c("ITC1179-Monyet", "Kokopo2")]




dat2 <- updog:::format_multidog(res)



