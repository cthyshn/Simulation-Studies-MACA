setwd(paste0("/bashscripts/diffmaf/", split))
  start <- cvsall$pos[i] - 125000
  stop <- cvsall$pos[i] + 125000
  
  file_name <- paste0("SuSiEx_", cvsall$id[i], ".sh")
  file_conn <- file(file_name, "w")
  
  writeLines(c(
    "#!/bin/bash -l",
    "#$ -P t2d-intern"), con = file_conn)
  writeLines(paste0("#$ -N SuSiEx_", gsub(":", "_", cvsall$id[i])), con = file_conn) 
  writeLines(c("#$ -j y",
    "#$ -m bae",
    "#$ -M cathy.shen@mail.mcgill.ca"), con = file_conn)
  
  dir_string <- sprintf("dir=\"/diffmaf/results/%s/%s/coloc+SuSiEx\"",  split, cvsall$id[i])
  writeLines(dir_string, con = file_conn)
  writeLines("cd \"$dir\"", con = file_conn)
  
  
  writeLines(sprintf("
for i in {1..100}
do
    out_name=\"gwas${i}.SuSiEx.EUR.AFR.output.cs95\"
    
    /SuSiEx/bin/SuSiEx --sst_file=EUR_gwas${i}.txt,AFR_gwas${i}.txt \\
    --n_gwas=%d,%d \\
    --ld_file=EUR,AFR \\
    --ref_file=EUR,AFR \\
    --plink=/SuSiEx/utilities/plink \\
    --out_dir=./ \\
    --out_name=$out_name \\
    --level=0.95 \\
    --pval_thresh=1 \\
    --maf=0.005 \\
    --chr=2 --bp=%d,%d \\
    --snp_col=2,2 \\
    --chr_col=1,1 \\
    --bp_col=3,3 \\
    --a1_col=4,4 \\
    --a2_col=5,5 \\
    --se_col=7,7 \\
    --pval_col=9,9 \\
    --eff_col=6,6 \\
    --min_purity=0.5 \\
    --n_sig=1 \\
    --keep-ambig=True \\
    --mult-step=True
done", eurn, afrn, start, stop), con = file_conn)
  
  writeLines(sprintf("
for i in {1..100}
do
    out_name=\"eqtl${i}.SuSiEx.EUR.AFR.output.cs95\"
    
    /SuSiEx/bin/SuSiEx --sst_file=EUR_eqtl${i}.txt,AFR_eqtl${i}.txt \\
    --n_gwas=%d,%d \\
    --ld_file=EUR,AFR \\
    --ref_file=EUR,AFR \\
    --plink=/SuSiEx/utilities/plink \\
    --out_dir=./ \\
    --out_name=$out_name \\
    --level=0.95 \\
    --pval_thresh=1 \\
    --maf=0.005 \\
    --chr=2 --bp=%d,%d \\
    --snp_col=2,2 \\
    --chr_col=1,1 \\
    --bp_col=3,3 \\
    --a1_col=4,4 \\
    --a2_col=5,5 \\
    --se_col=7,7 \\
    --pval_col=9,9 \\
    --eff_col=6,6 \\
    --min_purity=0.5 \\
    --n_sig=1 \\
    --keep-ambig=True \\
    --mult-step=True
done", eurn, afrn, start, stop), con = file_conn)
  
  
  close(file_conn)
  
  file_name <- paste0("msecaviar_", cvsall$id[i], ".sh")
  file_conn <- file(file_name, "w")
  

  
  writeLines(c(
    "#!/bin/bash -l",
    "#$ -P t2d-intern"), con = file_conn)
  writeLines(paste0("#$ -N msecaviar_", gsub(":", "_", cvsall$id[i])), con = file_conn) 
  writeLines(c("#$ -j y",
               "#$ -m bae",
               "#$ -M cathy.shen@mail.mcgill.ca"), con = file_conn)
  

  dir_string <- sprintf("dir=\"/diffmaf/results/%s/%s/MseCAVIAR\"", split, cvsall$id[i])
  writeLines(dir_string, con = file_conn)
  writeLines("cd \"$dir\"", con = file_conn)
  
  
  writeLines(sprintf("for i in {1..100}
    do
    ms_out_name=\"gwas${i}_eqtl${i}.msecaviar.EUR.AFR\"
    /eCAVIAR+MsCAVIAR/caviar/CAVIAR-C++/eCAVIAR -l ldfiles.txt -k ldfiles.txt -z gwas${i}.txt  -y eqtl${i}.txt -o $ms_out_name -r 0.95 -c 1 -n %d,%d -m %d,%d
    done", afrn, eurn, afrn, eurn), con = file_conn) 
  
  close(file_conn)
#}




