# Look at T2D by adiposity across ancestries

1. Create a cohort for cases and a cohort for controls on DNANexus
    1. In Cohort Browser, create a cohort that makes a filter for ICD10 code E11 (`Diagnoses - ICD10 includes E11 Non-insulin dependent diabetes mellitus`) and save it
    2. Create a comparison cohort as the inverse
2. Use Table Exporter to export the data. For each of the two cohorts
    - upload a file that has the field ids (see below)
    - Choose the cohort
    - add the field_id.txt file
    - Entity = `participant`
    - Set Output prefix


## Field IDs

```
eid
p23127_i0
p23099_i0
p21001_i0
p31
p21022
p21003_i0
```



