# Class 18: Pertussis and the CMI-PB project
Sofia Jaravata (A19160915)

- [Background](#background)
- [1. Investigating pertussis cases by
  year](#1-investigating-pertussis-cases-by-year)
  - [Q1.](#q1)
- [2. A tale of two vaccines (wP &
  aP)](#2-a-tale-of-two-vaccines-wp--ap)
  - [Q2.](#q2)
  - [Q3.](#q3)
- [3. Exploring CMI-PB data](#3-exploring-cmi-pb-data)
  - [The CMI-PB API returns JSON
    data](#the-cmi-pb-api-returns-json-data)
    - [Q4.](#q4)
    - [Q5.](#q5)
    - [Q6.](#q6)
  - [Side-Note: Working with dates](#side-note-working-with-dates)
    - [Q7.](#q7)
    - [Q8.](#q8)
    - [Q9.](#q9)
  - [Joining multiple tables](#joining-multiple-tables)
    - [Q9a.](#q9a)
    - [Q10.](#q10)
    - [Q11.](#q11)
    - [Q12.](#q12)
- [4. Examine IgG Ab titer levels](#4-examine-igg-ab-titer-levels)
  - [Q13.](#q13)
  - [Q14.](#q14)
  - [Q15.](#q15)
  - [Q16.](#q16)
  - [Q17.](#q17)
  - [Q18.](#q18)
- [5. Obtaining CMI-PB RNASeq data](#5-obtaining-cmi-pb-rnaseq-data)
  - [Q19.](#q19)
  - [Q20.:](#q20)
  - [Q21.](#q21)

# Background

Pertussis (more commonly known as **whooping cough**) is a highly
contagious respiratory disease caused by the bacterium Bordetella
pertussis (Figure 1). People of all ages can be infected leading to
violent coughing fits followed by a characteristic high-pitched “whoop”
like intake of breath. Children have the highest risk for severe
complications and death. Recent estimates from the WHO indicate that ~16
million cases and 200,000 infant deaths are due to pertussis annually
(Black et al. 2010).

# 1. Investigating pertussis cases by year

### Q1.

> With the help of the R “addin” package datapasta assign the CDC
> pertussis case number data to a data frame called `cdc` and use
> **ggplot** to make a plot of cases numbers over time.

``` r
library(datapasta)
cdc <- read.csv("U.S. Reported Pertussis Cases_ 1922 - 2025.csv")
cdc
```

        Year Number.of.Reported.Pertussis.Cases Data.Status
    1   1922                            107,473   Finalized
    2   1923                            164,191   Finalized
    3   1924                            165,418   Finalized
    4   1925                            152,003   Finalized
    5   1926                            202,210   Finalized
    6   1927                            181,411   Finalized
    7   1928                            161,799   Finalized
    8   1929                            197,371   Finalized
    9   1930                            166,914   Finalized
    10  1931                            172,559   Finalized
    11  1932                            215,343   Finalized
    12  1933                            179,135   Finalized
    13  1934                            265,269   Finalized
    14  1935                            180,518   Finalized
    15  1936                            147,237   Finalized
    16  1937                            214,652   Finalized
    17  1938                            227,319   Finalized
    18  1939                            103,188   Finalized
    19  1940                            183,866   Finalized
    20  1941                            222,202   Finalized
    21  1942                            191,383   Finalized
    22  1943                            191,890   Finalized
    23  1944                            109,873   Finalized
    24  1945                            133,792   Finalized
    25  1946                            109,860   Finalized
    26  1947                            156,517   Finalized
    27  1948                             74,715   Finalized
    28  1949                             69,479   Finalized
    29  1950                            120,718   Finalized
    30  1951                             68,687   Finalized
    31  1952                             45,030   Finalized
    32  1953                             37,129   Finalized
    33  1954                             60,886   Finalized
    34  1955                             62,786   Finalized
    35  1956                             31,732   Finalized
    36  1957                             28,295   Finalized
    37  1958                             32,148   Finalized
    38  1959                             40,005   Finalized
    39  1960                             14,809   Finalized
    40  1961                             11,468   Finalized
    41  1962                             17,749   Finalized
    42  1963                             17,135   Finalized
    43  1964                             13,005   Finalized
    44  1965                              6,799   Finalized
    45  1966                              7,717   Finalized
    46  1967                              9,718   Finalized
    47  1968                              4,810   Finalized
    48  1969                              3,285   Finalized
    49  1970                              4,249   Finalized
    50  1971                              3,036   Finalized
    51  1972                              3,287   Finalized
    52  1973                              1,759   Finalized
    53  1974                              2,402   Finalized
    54  1975                              1,738   Finalized
    55  1976                              1,010   Finalized
    56  1977                              2,177   Finalized
    57  1978                              2,063   Finalized
    58  1979                              1,623   Finalized
    59  1980                              1,730   Finalized
    60  1981                              1,248   Finalized
    61  1982                              1,895   Finalized
    62  1983                              2,463   Finalized
    63  1984                              2,276   Finalized
    64  1985                              3,589   Finalized
    65  1986                              4,195   Finalized
    66  1987                              2,823   Finalized
    67  1988                              3,450   Finalized
    68  1989                              4,157   Finalized
    69  1990                              4,570   Finalized
    70  1991                              2,719   Finalized
    71  1992                              4,083   Finalized
    72  1993                              6,586   Finalized
    73  1994                              4,617   Finalized
    74  1995                              5,137   Finalized
    75  1996                              7,796   Finalized
    76  1997                              6,564   Finalized
    77  1998                              7,405   Finalized
    78  1999                              7,298   Finalized
    79  2000                              7,867   Finalized
    80  2001                              7,580   Finalized
    81  2002                              9,771   Finalized
    82  2003                             11,647   Finalized
    83  2004                             25,827   Finalized
    84  2005                             25,616   Finalized
    85  2006                             15,632   Finalized
    86  2007                             10,454   Finalized
    87  2008                             13,278   Finalized
    88  2009                             16,858   Finalized
    89  2010                             27,550   Finalized
    90  2011                             18,719   Finalized
    91  2012                             48,277   Finalized
    92  2013                             28,639   Finalized
    93  2014                             32,971   Finalized
    94  2015                             20,762   Finalized
    95  2016                             17,972   Finalized
    96  2017                             18,975   Finalized
    97  2018                             15,609   Finalized
    98  2019                             18,617   Finalized
    99  2020                              6,124   Finalized
    100 2021                              2,116   Finalized
    101 2022                              3,044   Finalized
    102 2023                              7,063   Finalized
    103 2024                             43,321 Provisional
    104 2025                             28,783 Provisional

``` r
cdc$Number.of.Reported.Pertussis.Cases <- as.numeric(gsub("," , "" , cdc$Number.of.Reported.Pertussis.Cases))
```

``` r
library(ggplot2)
```

    Warning: package 'ggplot2' was built under R version 4.5.2

``` r
ggplot(cdc) +
  aes(Year, Number.of.Reported.Pertussis.Cases) +
  geom_point() +
  geom_line() +
  labs(x = "Year", 
       y = "Number of cases", 
       title = "Pertussis Cases by Year (1922- )")
```

![](class18_files/figure-commonmark/unnamed-chunk-3-1.png)

# 2. A tale of two vaccines (wP & aP)

### Q2.

> Using the ggplot geom_vline() function add lines to your previous plot
> for the 1946 introduction of the wP vaccine and the 1996 switch to aP
> vaccine (see example in the hint below). What do you notice?

``` r
ggplot(cdc) +
  aes(Year, Number.of.Reported.Pertussis.Cases) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept=c(1946, 1996),
             linetype = 2, 
             col = c("blue", "red")
             ) + 
  labs(x = "Year", y = "Number of cases", title = "Pertussis Cases by Year (1922- )")
```

![](class18_files/figure-commonmark/unnamed-chunk-4-1.png)

I notice that after the introduction of the wP vaccine (blue line), the
number of observed cases started to decrease.

### Q3.

> Describe what happened after the introduction of the aP vaccine? Do
> you have a possible explanation for the observed trend?

After the introduction of the aP vaccine in 1996, the number of observed
cases started increasing slowly. A possible explanation could be
PCR-based testing turning more sensitive, vaccination hesitancy,
bacterial evolution, or waning immunity in adolescents with the newer aP
vaccine in comparison to the older wP vaccine.

# 3. Exploring CMI-PB data

## The CMI-PB API returns JSON data

``` r
# Allows us to read, write and process JSON data
library(jsonlite)
```

``` r
subject <- read_json("https://www.cmi-pb.org/api/subject", simplifyVector = TRUE) 
head(subject, 3)
```

      subject_id infancy_vac biological_sex              ethnicity  race
    1          1          wP         Female Not Hispanic or Latino White
    2          2          wP         Female Not Hispanic or Latino White
    3          3          wP         Female                Unknown White
      year_of_birth date_of_boost      dataset
    1    1986-01-01    2016-09-12 2020_dataset
    2    1968-01-01    2019-01-28 2020_dataset
    3    1983-01-01    2016-10-10 2020_dataset

### Q4.

> How many aP and wP infancy vaccinated subjects are in the dataset?

``` r
table(subject$infancy_vac)
```


    aP wP 
    87 85 

There are 87 aP and 85 wP infancy vaccinated subjects in the dataset.

### Q5.

> How many Male and Female subjects/patients are in the dataset?

``` r
table(subject$biological_sex)
```


    Female   Male 
       112     60 

There are 112 Female patients and 60 Male patients.

### Q6.

> What is the breakdown of race and biological sex (e.g. number of Asian
> females, White males etc…)?

``` r
table(subject$race, subject$biological_sex)
```

                                               
                                                Female Male
      American Indian/Alaska Native                  0    1
      Asian                                         32   12
      Black or African American                      2    3
      More Than One Race                            15    4
      Native Hawaiian or Other Pacific Islander      1    1
      Unknown or Not Reported                       14    7
      White                                         48   32

Breakdown of race and biological sex displayed in table above.

## Side-Note: Working with dates

### Q7.

> Using this approach determine (i) the average age of wP individuals,
> (ii) the average age of aP individuals; and (iii) are they
> significantly different?

``` r
library(lubridate)
```


    Attaching package: 'lubridate'

    The following objects are masked from 'package:base':

        date, intersect, setdiff, union

``` r
# Use todays date to calculate age in days
subject$age <- today() - ymd(subject$year_of_birth)
```

``` r
library(dplyr)
```

    Warning: package 'dplyr' was built under R version 4.5.2


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
#Split by aP group, summarize aP ages
ap <- subject %>% filter(infancy_vac == "aP")
round( summary( time_length( ap$age, "years" ) ) )
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
         23      27      28      28      29      35 

``` r
# Split by wP group, summarize wP ages
wp <- subject %>% filter(infancy_vac == "wP")
round(summary(time_length( wp$age, "years" )))
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
         23      33      35      37      40      58 

The average age of wP individuals is 37 and 28 for aP individuals. Based
on the results, it does seem that they are significantly different as
the wP values are greater, and the Median, Mean, 3rd Qu., and Max wP
values different and much greater than the aP values.

### Q8.

> Determine the age of all individuals at time of boost?

``` r
int <- ymd(subject$date_of_boost) - ymd(subject$year_of_birth)
age_at_boost <- time_length(int, "year")
head(age_at_boost)
```

    [1] 30.69678 51.07461 33.77413 28.65982 25.65914 28.77481

### Q9.

> With the help of a faceted boxplot or histogram (see below), do you
> think these two groups are significantly different?

``` r
ggplot(subject) +
  aes(time_length(age, "year"),
      fill=as.factor(infancy_vac)) +
  geom_histogram(show.legend=FALSE) +
  facet_wrap(vars(infancy_vac), nrow=2) +
  xlab("Age in years")
```

    `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

![](class18_files/figure-commonmark/unnamed-chunk-15-1.png)

``` r
# Or use wilcox.test() 
x <- t.test(time_length( wp$age, "years" ),
       time_length( ap$age, "years" ))

x$p.value
```

    [1] 2.372101e-23

I think these two groups are significantly different as the p-value of
2.372101e-23 indicates that they are different. The aP group has more
cases focused in those who are younger with ages ~25-30, while the wP
group has more cases spread across all ages, mainly focused in the age
range of 30-45, with some cases found in ages 45+.

## Joining multiple tables

Read the specimen and ab_titer tables into R and store the data as
specimen and titer named data frames.

``` r
# Complete the API URLs...
specimen <- read_json("https://www.cmi-pb.org/api/specimen", simplifyVector = TRUE) 
titer <- read_json("https://www.cmi-pb.org/api/plasma_ab_titer", simplifyVector = TRUE) 
```

### Q9a.

> Complete the code to join specimen and subject tables to make a new
> merged data frame containing all specimen records along with their
> associated subject details:

``` r
meta <- inner_join(specimen, subject)
```

    Joining with `by = join_by(subject_id)`

``` r
dim(meta)
```

    [1] 1503   14

``` r
head(meta)
```

      specimen_id subject_id actual_day_relative_to_boost
    1           1          1                           -3
    2           2          1                            1
    3           3          1                            3
    4           4          1                            7
    5           5          1                           11
    6           6          1                           32
      planned_day_relative_to_boost specimen_type visit infancy_vac biological_sex
    1                             0         Blood     1          wP         Female
    2                             1         Blood     2          wP         Female
    3                             3         Blood     3          wP         Female
    4                             7         Blood     4          wP         Female
    5                            14         Blood     5          wP         Female
    6                            30         Blood     6          wP         Female
                   ethnicity  race year_of_birth date_of_boost      dataset
    1 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    2 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    3 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    4 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    5 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    6 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
             age
    1 14759 days
    2 14759 days
    3 14759 days
    4 14759 days
    5 14759 days
    6 14759 days

### Q10.

> Now using the same procedure join meta with titer data so we can
> further analyze this data in terms of time of visit aP/wP, male/female
> etc.

``` r
abdata <- inner_join(titer, meta)
```

    Joining with `by = join_by(specimen_id)`

``` r
dim(abdata)
```

    [1] 52576    21

### Q11.

> How many specimens (i.e. entries in abdata) do we have for each
> isotype?

``` r
table(abdata$isotype)
```


      IgE   IgG  IgG1  IgG2  IgG3  IgG4 
     6698  5389 10117 10124 10124 10124 

We have **6698** specimens for IgE, **5389** for IgG, **10117** for IgG1
and **10124** for IgG2, IgG3, IgG4.

### Q12.

> What are the different \$dataset values in abdata and what do you
> notice about the number of rows for the most “recent” dataset?

``` r
table(abdata$dataset)
```


    2020_dataset 2021_dataset 2022_dataset 2023_dataset 
           31520         8085         7301         5670 

The different \$dataset values in abdata are 31520 (2020 dataset), 8085
(2021 dataset), 7301 (2022 dataset), and 5670 (2023 dataset). The number
of rows for the most “recent” dataset, 5670 in the 2023 dataset, is much
lower than the other datasets, especially in comparison to the number of
rows in the 2020 dataset, 31520.

# 4. Examine IgG Ab titer levels

Now using our joined/merged/linked `abdata` dataset `filter(`) for IgG
`isotype`.

``` r
igg <- abdata %>% filter(isotype == "IgG")
head(igg)
```

      specimen_id isotype is_antigen_specific antigen        MFI MFI_normalised
    1           1     IgG                TRUE      PT   68.56614       3.736992
    2           1     IgG                TRUE     PRN  332.12718       2.602350
    3           1     IgG                TRUE     FHA 1887.12263      34.050956
    4          19     IgG                TRUE      PT   20.11607       1.096366
    5          19     IgG                TRUE     PRN  976.67419       7.652635
    6          19     IgG                TRUE     FHA   60.76626       1.096457
       unit lower_limit_of_detection subject_id actual_day_relative_to_boost
    1 IU/ML                 0.530000          1                           -3
    2 IU/ML                 6.205949          1                           -3
    3 IU/ML                 4.679535          1                           -3
    4 IU/ML                 0.530000          3                           -3
    5 IU/ML                 6.205949          3                           -3
    6 IU/ML                 4.679535          3                           -3
      planned_day_relative_to_boost specimen_type visit infancy_vac biological_sex
    1                             0         Blood     1          wP         Female
    2                             0         Blood     1          wP         Female
    3                             0         Blood     1          wP         Female
    4                             0         Blood     1          wP         Female
    5                             0         Blood     1          wP         Female
    6                             0         Blood     1          wP         Female
                   ethnicity  race year_of_birth date_of_boost      dataset
    1 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    2 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    3 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    4                Unknown White    1983-01-01    2016-10-10 2020_dataset
    5                Unknown White    1983-01-01    2016-10-10 2020_dataset
    6                Unknown White    1983-01-01    2016-10-10 2020_dataset
             age
    1 14759 days
    2 14759 days
    3 14759 days
    4 15855 days
    5 15855 days
    6 15855 days

### Q13.

> Complete the following code to make a summary boxplot of Ab titer
> levels (MFI) for all antigens:

``` r
ggplot(igg) +
  aes(MFI_normalised, antigen, col = infancy_vac) +
  geom_boxplot(show.legend = FALSE) + 
    xlim(0,75) +
  facet_wrap(vars(visit), nrow=2)
```

    Warning: Removed 5 rows containing non-finite outside the scale range
    (`stat_boxplot()`).

![](class18_files/figure-commonmark/unnamed-chunk-23-1.png)

### Q14.

> What antigens show differences in the level of IgG antibody titers
> recognizing them over time? Why these and not others?

The FHA and FIM2/3 antigens show differences in the level of IgG
antibody titers recognizing them, as well as little changes in PT and
PRN, but not as much change in TT, DT, and OVA. This is because these
antigens, FHA, FIM2/3, PT, and PRN, are Bordetella pertussis antigens,
and they exhibit significantly elevated antibody titers, demonstrating
an effective immune response directed at pertussis-specific antigens.
OVA is used as a negative control antigen and thus not targeted by the
vaccine and TT and DT are more stable comparison antigens. If the immune
response is changing after pertussi exposure or vaccination, the
antigens showing differences would most likely show it.

### Q15.

> Filter to pull out only two specific antigens for analysis and create
> a boxplot for each. You can chose any you like. Below I picked a
> “control” antigen (“OVA”, that is not in our vaccines) and a clear
> antigen of interest (“PT”, Pertussis Toxin, one of the key virulence
> factors produced by the bacterium B. pertussis).

``` r
filter(igg, antigen=="OVA") %>%
  ggplot() +
  aes(MFI_normalised, col=infancy_vac) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(vars(visit)) +
  theme_bw()+
  labs(title = "OVA antigen levels per visit (aP red, wP teal)")
```

![](class18_files/figure-commonmark/unnamed-chunk-24-1.png)

***and the same for antigen==“PT”***

``` r
filter(igg, antigen=="PT") %>%
  ggplot() +
  aes(MFI_normalised, col=infancy_vac) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(vars(visit)) +
  theme_bw() +
  labs(title = "PT antigen levels per visit (aP red, wP teal)")
```

![](class18_files/figure-commonmark/unnamed-chunk-25-1.png)

### Q16.

> What do you notice about these two antigens time courses and the PT
> data in particular?

As time passes, OVA levels stay pretty consistently low among all visits
as OVA is the negative control antigen. Over time, PT antibody levels
are low overall during visits 1-3, start to increase in visits 4-7, and
then decreases to baseline for visits 8-12, indicating that the PT
response is variable. This comparison of OVA vs. PT shows that
antigen-specific immunity grows against the vaccine.

### Q17.

> Do you see any clear difference in aP vs. wP responses?

There are not many clear differences in aP vs. wP responses, they seem
pretty similar. Both groups increase in the PT antibody titers, but
eventually decrease. Thus, both vaccine types can result in immune
responses.

**Lets finish this section by looking at the 2021 dataset IgG PT antigen
levels time-course:**

``` r
abdata.21 <- abdata %>% filter(dataset == "2021_dataset")

abdata.21 %>% 
  filter(isotype == "IgG",  antigen == "PT") %>%
  ggplot() +
    aes(x=planned_day_relative_to_boost,
        y=MFI_normalised,
        col=infancy_vac,
        group=subject_id) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_vline(xintercept=14, linetype="dashed") +
  labs(title="2021 dataset IgG PT",
       subtitle = "Dashed lines indicate day 0 (pre-boost) and 14 (apparent peak levels)")
```

![](class18_files/figure-commonmark/unnamed-chunk-26-1.png)

``` r
abdata.20 <- abdata %>% filter(dataset == "2020_dataset")

abdata.20 %>% 
  filter(isotype == "IgG",  antigen == "PT") %>%
  ggplot() +
    aes(x=planned_day_relative_to_boost,
        y=MFI_normalised,
        col=infancy_vac,
        group=subject_id) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_vline(xintercept=14, linetype="dashed") +
  labs(title="2020 dataset IgG PT",
       subtitle = "Dashed lines indicate day 0 (pre-boost) and 14 (apparent peak levels)")
```

![](class18_files/figure-commonmark/unnamed-chunk-27-1.png)

### Q18.

> Does this trend look similar for the 2020 dataset?

Yes, this trend looks pretty similar for the 2020 dataset. The IgG PT
levels are peak around day 14 before their decline.

# 5. Obtaining CMI-PB RNASeq data

``` r
url <- "https://www.cmi-pb.org/api/v2/rnaseq?versioned_ensembl_gene_id=eq.ENSG00000211896.7"

rna <- read_json(url, simplifyVector = TRUE) 
```

``` r
#meta <- inner_join(specimen, subject)
ssrna <- inner_join(rna, meta)
```

    Joining with `by = join_by(specimen_id)`

### Q19.

> Make a plot of the time course of gene expression for IGHG1 gene
> (i.e. a plot of visit vs. tpm).

``` r
ggplot(ssrna) +
  aes(visit, tpm, group=subject_id) +
  geom_point() +
  geom_line(alpha=0.2)
```

![](class18_files/figure-commonmark/unnamed-chunk-30-1.png)

### Q20.:

> What do you notice about the expression of this gene (i.e. when is it
> at it’s maximum level)?

I notice that the expression of this gene is at its maximum level right
after booster vaccination. The IGHG1 expression is mostly low at many
visits, but then spikes early with a maximum level at visit 4 and then
randomly later on at visit 8.

### Q21.

> Does this pattern in time match the trend of antibody titer data? If
> not, why not?

No, this pattern in time does not match the trend of antibody titer data
because the IGHG1 expression level peaks at visit 4, while the antibody
titer data displays variable data and peaks later on. IGHG1 expression
measures mRNA inside cells and antibody titer measures secreted antibody
protein in blood. IGHG1 expression levels peak earlier and antibody
titers peak later because the mRNA must first be translated into
proteins before those antibodies are secreted into the bloodstream.

*We can dig deeper and color and/or facet by infancy_vac status:*

``` r
ggplot(ssrna) +
  aes(tpm, col=infancy_vac) +
  geom_boxplot() +
  facet_wrap(vars(visit))
```

![](class18_files/figure-commonmark/unnamed-chunk-31-1.png)

There is however no obvious wP vs. aP differences here even if we focus
in on a particular visit:

``` r
ssrna %>%  
  filter(visit==4) %>% 
  ggplot() +
    aes(tpm, col=infancy_vac) + geom_density() + 
    geom_rug() 
```

![](class18_files/figure-commonmark/unnamed-chunk-32-1.png)
