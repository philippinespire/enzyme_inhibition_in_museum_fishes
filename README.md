# enzyme_inhibition_in_museum_fishes

This repo contains the scripts used in the study:
  
  ***"Effectiveness of Reduced Representation Sequencing on Century-Old Ethanol-Preserved Museum Fishes"*** 

which compares error rates and read yeild from RADseq libraries constructed from contemporary and ethanol-preserved marine fishes collected over 100 years ago.  Compared to contemporary samples, results from historical specimens illustrated elevated read loss and error rates in both the synthetic adapter sequence and the last two positions of the Sbf1 restriction site (natural fish sequence) indicating that the specificity of both the DNA polymerase and restriction enzyme, respectively, was impaired. We hypothesize that the observed enzyme inhibition is caused by an unknow agent in the preservative of the museum specimens. Identifying this agent will require additional studies.

Overall, sequencing of fishes preserved and stored in EtOH for >100 years is possible, but all else being equal, it can result in more sequence substitution errors, unintended loci, and decreased depth of coverage due to altered enzymatic activity during library preparation when compared to contemporary samples. Consequently, up to 24% more DNA per museum specimen needed to be sequenced to achieve comparable results to contemporary specimens in the current study.

---
### Citation

Effectiveness of Reduced Representation Sequencing on Century-Old Ethanol-Preserved Museum Fishes (2024) Molecular Ecology Resources

---
## Scripts

<details><summary>ChiSQTest.R</summary>
<p>

[ChiSQTest.R](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/ChiSQTest.R)

</p>
</details>

<details><summary>countReads.sbatch</summary>
<p>

</p>
</details>
<details><summary>ErrorRateCalc.bash</summary>
<p>
[ErrorRateCalc.bash](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/ErrorRateCalc.bash)

</p>
</details>
<details><summary>ErrorRateCodeFunction</summary>
<p>
[ErrorRateCodeFunction](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/ErrorRateCodeFunction)

</p>
</details>
<details><summary>mapDamage.sh</summary>
<p>

[mapDamage.sh](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/mapDamage.sh)

</p>
</details>
<details><summary>PIRE_Stats_forR.R</summary>
<p>

[PIRE_Stats_forR.R](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/PIRE_Stats_forR.R)
</p>
</details>
<details><summary>processSequenceCounts.R</summary>
<p>

[processSequenceCounts.R](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/processSequenceCounts.R)

</p>
</details>
<details><summary>SequenceCounts.bash

[SequenceCounts.bash](https://github.com/philippinespire/enzyme_inhibition_in_museum_fishes/blob/main/SequenceCounts.bash)
</p>
</details>
