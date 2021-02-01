# Swiss data

Individual data on cases of SARS-CoV-2 infection are released regularly by the Federal Office of Public Health (FOPH, https://www.bag.admin.ch/bag/en/home.html) to Swiss scientists (not to be shared with the public). Data on all-causes mortality are released by the Federal Statistical Office (FSO, https://www.bfs.admin.ch/bfs/en/home.html) via the same channel. 

## Dictionary

### Individual data: `indiv_data`

Anonymized data about every case of SARS-CoV-2 infection reported in Switzerland (n=30,691 on 25 May 2020). 

- `id`: unique identifier
- `ktn`: canton in 2 letters (26 cantons plus Liechtenstein)
- `sex`: male or female
- `age_class`: nine age classes from 0-9 to 80+
- `onset_dt`: date of disease onset (available only for a subset)
- `test_dt`: date of testing
- `report_dt`: date of reporting the positive test to the FOPH
- `hospit`: hospitalisation (0=no, 1=yes)
- `hospit_dt`: date of hospitalisation
- `death`: death (0=no, 1=yes)
- `death_dt`: date of death
- `onset_to_test`: delay in days from onset to test
- `(...)`: all delays between events
- `ktn_id`: numeric (1:27)
- `canton`: full canton name (26 cantons plus Liechtenstein)
- `region`: NUTS-2 region (7 regions, Liechtenstein is put with Ostschweiz)
- `region_id`: numeric (1:7)

### Data aggregated by date nationally: `agg_data`

- `date`: from 1 January 2020 to 25 May 2020
- `onset_dt`: number of cases by date of disease onset (for the subset with this information available)
- `test_dt`: number of cases by date of testing
- `report_dt`: number of cases by date of report
- `hospit_dt`: number of hospitalisation admissions by date
- `death_dt`: number of deaths by date
- `death_dt_hospit`: number of deaths with hospitalisation by date
- `death_dt_nohospit`: number of deaths without hospitalisation by date
- `total_tested`: number of tests performed by date

### Data aggregated by date and canton: `agg_data_ktn`

- `ktn_id`: numeric (1:27)
- `date`: from 1 January 2020 to 25 May 2020
- `onset_dt`: number of cases by date of disease onset (for the subset with this information available)
- `test_dt`: number of cases by date of testing
- `report_dt`: number of cases by date of report
- `hospit_dt`: number of hospitalisation admissions by date
- `death_dt`: number of deaths by date
- `death_dt_hospit`: number of deaths with hospitalisation by date
- `death_dt_nohospit`: number of deaths without hospitalisation by date
- `total_tested`: NOT AVAILABLE YET
- `ktn`: canton in 2 letters (26 cantons plus Liechtenstein)
- `canton`: full canton name (26 cantons plus Liechtenstein)
- `region`: NUTS-2 region (7 regions, Liechtenstein is put with Ostschweiz)
- `region_id`: numeric (1:7)
- `pop`: canton population (2019)


### Data aggregated by date and region: `agg_data_region`

- `region_id`: numeric (1:7)
- `date`: from 1 January 2020 to 25 May 2020
- `onset_dt`: number of cases by date of disease onset (for the subset with this information available)
- `test_dt`: number of cases by date of testing
- `report_dt`: number of cases by date of report
- `hospit_dt`: number of hospitalisation admissions by date
- `death_dt`: number of deaths by date
- `death_dt_hospit`: number of deaths with hospitalisation by date
- `death_dt_nohospit`: number of deaths without hospitalisation by date
- `total_tested`: NOT AVAILABLE YET
- `region`: NUTS-2 region (7 regions, Liechtenstein is put with Ostschweiz)
- `pop_region`: region population (2019)

### Individual data on all-causes mortality in Switzerland: `all_causes_mortality`

- `death_dt`: from 1 January 2020 to 22 May 2020
- `age_class`: nine age classes from 0-9 to 80+
- `sex`: male or female
- `region_id`: numeric (1:5)
- `region`: NUTS-2 region (7 regions, Liechtenstein is put with Ostschweiz)
- `pop_region`: region population (2019)

### Meta-data: `ktn_key` and `region_key`


