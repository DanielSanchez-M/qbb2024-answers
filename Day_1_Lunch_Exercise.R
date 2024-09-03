#import tidyverse packages
library("tidyverse")

#3. Load File and assign to df
df <- read_tsv("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

#3. Create Subject column
df <- df %>%
  mutate( SUBJECT=str_extract( SAMPID, "[^-]+-[^-]+" ), .before=1)

#4. Use group_by, summarize, and arrange to see which two SUBJECTS have most samples
df %>%
  group_by(SUBJECT) %>%
  summarize(n_of_samples=n()) %>%
  arrange(desc(n_of_samples))
##most samples = K-562 & GTEX-NPJ8
df %>%
  group_by(SUBJECT) %>%
  summarize(n_of_samples=n()) %>%
  arrange(n_of_samples)
##least samples = GTEX-1JMI6 & GTEX-1PAR6

#5. Which two SMTSDs have most samples? Least?
df %>%
  group_by(SMTSD) %>%
  summarize(n_of_samples=n()) %>%
  arrange(desc(n_of_samples))
##most samples = Whole Blood, Muscle-Skeletal
##Why? = Blood is easier to collect, and theres an abundance of them
df %>%
  group_by(SMTSD) %>%
  summarize(n_of_samples=n()) %>%
  arrange(n_of_samples)
##least samples = Kidney-Medula, Cervix-Ectocervix/Fallopian Tube
##Why? = Minimal # of tissues, poor integrity

#6. For GTEX-NPJ8: Filter and Save new Object.
df_npj8 <- df %>%
  filter(SUBJECT == "GTEX-NPJ8")
 
df_npj8 %>%
  group_by(SMTSD) %>%
  summarize(counts=n()) %>%
  arrange(-counts)
##Most Samples = Whole Blood - 9
## Difference = They used different sequencing techniques

#7. SMATSSCR
df %>%
  filter(!is.na(SMATSSCR)) %>%
  group_by(SMATSSCR) %>%
  summarize(n())
##A majority of the samples have a mean autolysis score of 1 and fewer samples have an autolysis score greater than or equal to 2. 0 = 3554, 1= 10,410, 2 = 1582, 3 = 193
## This report can be presented as a histogram or a bar graph