# DFCI G2C Database

Copyright (c) 2023-Present, [Ryan L. Collins](mailto:Ryan_Collins@dfci.harvard.edu) and the [Collins Genomics Lab](https://labs.dana-farber.org/collins-genomics/) and [Van Allen Lab](https://labs.dana-farber.org/vanallenlab/) at the Dana-Farber Cancer Institute.  

Distributed under terms of the [GNU GPL v2.0 License](/LICENSE) (see `LICENSE`).  


---  

This subdirectly contains code to curate participant- and sample-level metadata for all cohorts included in G2C.  

Given the highly heterogeneous nature of the cohorts included in G2C, we have mapped all available phenotype information onto a simplified data dictionary as follows:  

1. `reported_sex` : reported sex of participant, if known. Eligible values are `male`, `female`, `other`, or `NA` if not reported. Note that this is not necessarily consistent with inferred sex from WGS ploidy.  

2. `reported_race_or_ethnicity` : participant race and/or ethnicity, if either are reported. Note that this will frequently differ from WGS-inferred genetic ancestry. In some cases (e.g., certain TOPMed cohorts, 1000 Genomes Project), genetic ancestry inference had been previously performed, in which case these values superceded self-reported race and/or ethnicity. `NA` if not reported.    

3. `age` : individual age in years, if known. Ideally this reflects age at first primary cancer diagnosis, if reported, otherwise this value should reflect age at study enrollment or consent. Otherwise, this is assumed to correspond to the time of sample acquisition for sequencing (e.g., time at blood draw used for WGS). `NA` if unknown.    

4. `birth_year` : year of birth; `NA` if unknown.  

5. `vital_status` : one-hot indicator for participant vital status, if known. Eligible values are `0` if individual is known to have deceased since recruitment or diagnosis, `1` if individual is known to be alive at last contact, or `NA` if no follow-up was conducted, vital status is unknown, or otherwise not reported.  

6. `age_at_last_contact` : age in years at time of last contact or death. If no follow-up was conducted or reported, this value will match `age`.  

7. `years_to_last_contact` : time elapsed from `age` to most recent contact or death, if reported. Ideally, this value should reflect time from initial cancer diagnosis, otherwise it should reflect follow-up time from initial participant recruitment into the study. Should be interpreted in conjunction with `years_left_ceonsored` to avoid immortal time bias. `NA` if no follow-up was conducted or value is otherwise unknown or not reported.  

8. `years_left_censored` : unreportable/immortal time from `age` to `age_at_last_contact`, usually due to study enrollment or cancer diagnosis predating sample acquisition for sequencing. Should be evaluated jointly with `years_to_last_contact`. `NA` if no follow-up was conducted or value is otherwise unknown or not reported.  

9. `height` : participant height in centimeters, if reported, presumably at time of recruitment. `NA` if not reported.  

10. `weight` : participant height in kilograms, if reported, presumably at time of recruitment. `NA` if not reported.  

11. `bmi` : participant BMI, if reported or if both `height` and `weight` are available, presumably at time of recruitment. `NA` if not reported or not able to be computed.  

12. `cancer` : primary cancer diagnos(es) for each patient. Metastatic/secondary sites are excluded from consideration in these fields where clinical records are sufficiently detailed to allow disambiguation. This field can have one or more values delimited with a semicolon; note that in most cases there will be a single reported primary cancer diagnosis for a patient. Diagnoses were collapsed into the following 16 domains: `prostate`, `breast`, `lung`, `colorectal`, `melanoma`, `uterus`, `kidney`, `bladder`, `oral_cavity`, `ovary`, `cns`, `pancreas`, `esophagus`, `liver`, `stomach`, `other`. Note that `other` is a catch-all category for either (i) a known cancer diagnosis that does not map onto any of the 15 primary cancer domains or (ii) the individual is known to have a cancer diagnosis but the nature of the cancer is not reported. Cancer-free controls have `control` in this field; note that some phenotyping error is expected and permitted in this definition, but by-and-large we expect `control` samples to have rates of cancer at or below the overall rates in the general population. Samples with inadequate information to determine their (likely/plausible) cancer status will have `unknown` in this field, but wherever possible another label (e.g., `control` or `other`) is preferred. `NA` is not permitted.  

13. `stage` : AJCC disease stage for first primary cancer, if reported. Ideally recorded at time of diagnosis. Pathology staging preferred, followed by clinical staging, followed by pathology-reported T category. Eligible values are `0`, `I`, `II`, `III`, `IV`, `unknown` if not reported, and `NA` for cancer-free controls.  

14. `metastatic` : one-hot indicator if patient is known to have developed metastatic disease. Eligible values are `0`, `1`, `unknown` if not reported, and `NA` for cancer-free controls. Note that this value may disagree with `stage` if the patient was staged at diagnosis then subsequently developed metastatic disease.  

15. `grade` : Grade of primary tumor, if reported. Ideally reported at time of diagnosis. Eligible values are `1`, `2`, `3`, `4`, `unknown` if not reported, and `NA` for cancer-free controls. Borderline grade tumors (`GB`) were treated as grade `1`.  

16. `smoking_history` : one-hot indicator for patient smoking history, if known. Eligible values include `0` for reported never smoker, `1` for positive history of smoking, and `NA` for unreported or missing data.  

17. `cancer_icd10` : semicolon-delimited ICD-10 codes for this patient's cancer diagnoses, if reported. `NA` if missing and for cancer-free controls.  

18. `original_dx` : free text field to record original specific cancer diagnosis, which is not standardized across cohorts. `NA` for missing values and cancer-free controls.  

19. `wgs_tissue` : DNA source (i.e., tissue or biospecimen) used for germline whole-genome sequencing. `unknown` if not reported; `NA` not allowed.  

---  

Note that access to individual level phenotype data records are protected by agreements and consents from various entities (e.g., dbGaP, NIH All of Us, Hartwig Medical Foundation), so no blanket access can be provided. The above data descriptors and scripts contained in this subdirectory are sufficient to recreate the phenotype data used in all analyses, however.  

