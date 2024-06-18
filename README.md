# enzyme_inhibition_in_museum_fishes

This repo contains the scripts used in the study:
  
  ***"Effectiveness of Reduced Representation Sequencing on Century-Old Ethanol-Preserved Museum Fishes"*** 

which compares error rates and read yeild from RADseq libraries constructed from contemporary and ethanol-preserved marine fishes collected over 100 years ago.  Compared to contemporary samples, results from historical specimens illustrated elevated read loss and error rates in both the synthetic adapter sequence and the last two positions of the Sbf1 restriction site (natural fish sequence) indicating that the specificity of both the DNA polymerase and restriction enzyme, respectively, was impaired. We hypothesize that the observed enzyme inhibition is caused by an unknow agent in the preservative of the museum specimens. Identifying this agent will require additional studies.

Overall, sequencing of fishes preserved and stored in EtOH for >100 years is possible, but all else being equal, it can result in more sequence substitution errors, unintended loci, and decreased depth of coverage due to altered enzymatic activity during library preparation when compared to contemporary samples. Consequently, up to 24% more DNA per museum specimen needed to be sequenced to achieve comparable results to contemporary specimens in the current study.

---
### Citation

Effectiveness of Reduced Representation Sequencing on Century-Old Ethanol-Preserved Museum Fishes (2024) Molecular Ecology Resources

___ 
### Funding

This repo and forementioned publication are products of the collaborative [Philippines PIRE Project](https://sites.wp.odu.edu/PIRE/), [NSF Award #1743711](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1743711).

---
## Scripts

<details><summary>ChiSQTest.R</summary>
<p>

[ChiSQTest.R](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/ChiSQTest.R) performs data cleaning, transformation, and conducts chi-square tests to assess goodness-of-fit on sequencing error data from the two different fish species and time points. Following this, a Bayesian multinomial regression model is used to analyze error rates, fitting the model to the data and plotting the results. Finally, post-hoc tests are used to compare error rates between species and time points, providing a comprehensive statistical analysis of sequencing errors across different conditions.

</p>
</details>

<details><summary>countReads.sbatch</summary>
<p>

[countReads.sbatch](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/countReads.sbatch) provides code for counting reads from compressed FASTQ sequence files in parallel

</p>
</details>

<details><summary>ErrorRateCalc.bash</summary>
<p>

[ErrorRateCalc.bash](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/ErrorRateCalc.bash) calculates error rates from RAD sequencing data and uses the agrep command to identify specific patterns within the sequences and generate error-related data. 

</p>
</details>

<details><summary>ErrorRateCodeFunction</summary>
<p>

The [ErrorRateCodeFunction](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/ErrorRateCodeFunction) processes FASTQ files to calculate error rates associated with specific barcodes and reads. This function calculates the number of deletions and insertions at various positions within the reads using the command agrep, and computes the proportions of these outcomes relative to the total number of sequences. Additionally, the function calculates the proportion of substitutions at specific positions where there are no deletions. 

</p>
</details>

<details><summary>mapDamage.sh</summary>
<p>

[mapDamage.sh](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/mapDamage.sh) is a simple bash script to run mapDamage in BAM files in parallel

</p>
</details>

<details><summary>PIRE_Stats_forR.R</summary>
<p>

[PIRE_Stats_forR.R](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/PIRE_Stats_forR.R) processes, analyzes and plots sequencing error data. It calculates sequencing errors by comparing the base sequences, categorizes them into groups (barcode, ligation, and naturals sites), and computes error and success rates. The data is then summarized and plotted to visualize mean error rates and the distribution of insertions and deletions (Indels). The script also builds statistical models to estimate error rates, performs pairwise comparisons, and generates plots to compare observed means with model estimates. This script finally performs an analysis to predict error rates at position 15 and 16 using various predictors, and assesses the fit of these models.

</p>
</details>

<details><summary>processSequenceCounts.R</summary>
<p>

[processSequenceCounts.R](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/processSequenceCounts.R) processes and analyzes sequence yield data across species and treatments. The script generates bar plots of read proportions and counts, and conducts linear and beta regression analyses to explore relationships between species, collection time, and sequence counts. Finally, it runs mixed-effects models to test for differences in read proportions and counts.

</p>
</details>

<details><summary>SequenceCounts.bash</summary>
</p>
  
[SequenceCounts.bash](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/SequenceCounts.bash) counts read numbers in the input FASTQ files and perform subsequent analyses to categorize reads based on the number of mismatches from the expected sequences (2Del to 8Ins)

</p>
</details>
