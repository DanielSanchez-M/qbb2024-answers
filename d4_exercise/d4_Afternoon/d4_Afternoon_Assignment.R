
library(tidyverse)
library(broom)

#1.1 First, load aau1043_dnm.csv into a dataframe
dnm <- read_csv( file = "~/qbb2024-answers/d4_exercise/d4_Afternoon/aau1043_dnm.csv")

#1.2 tabulate the number of paternally and maternally inherited DNMs in each proband
dnm_summary <- dnm %>%
  group_by(Proband_id) %>%
  summarize(n_paternal_dnm = sum(Phase_combined == "father", na.rm = TRUE),
            n_maternal_dnm = sum(Phase_combined == "mother", na.rm = TRUE))
print(dnm_summary)

#1.3 load the data from aau1043_parental_age.csv into a new dataframe
ages <- read_csv( file = "~/qbb2024-answers/d4_exercise/d4_Afternoon/aau1043_parental_age.csv")

#1.4 Use the left_join() function to combine your dataframe from step 2 with the dataframe you just created in step 3 based on the shared column Proband_id
dnm_by_parental_age <- left_join(dnm_summary, ages, by = "Proband_id")

#2.1 exploring if there’s a relationship between the number of DNMs and parental age.
#   1. the count of maternal de novo mutations vs. maternal age
ggplot(data = dnm_by_parental_age,
       mapping = aes(
         x = Mother_age,
         y = n_maternal_dnm)) +
  geom_point(size = 0.75) + 
  geom_smooth(method = "lm", size = .75) +
  labs(x = "Maternal Age", y = "Number of Maternal de novo mutations")
#   2. the count of paternal de novo mutations vs. paternal age
ggplot(data = dnm_by_parental_age,
       mapping = aes(
         x = Father_age,
         y = n_paternal_dnm)) +
  geom_point(size = 0.75) + 
  geom_smooth(method = "lm", size = .75) +
  labs(x = "Paternal Age", y = "Number of Paternal de novo mutations")

#2.2 Fit a linear regression model to the data using the lm() function for Maternal.
lm( data = dnm_by_parental_age,
    formula = n_maternal_dnm ~ 1 + Mother_age) %>%
  summary()
#   1.1 What is the “size” of this relationship?
##    The size of this relationship is 0.37757
#   1.2 In your own words, what does this mean?
##    The coefficient of 0.37757 indicates that as the mother's age increases by 1 year, the number of maternal DNMs is expected to increase by around 0.38.
#   1.3 Does this match what you observed in your plots in step 2.1?
##    Yes, we can observe a positive trend line in the data

#   2.1 Is this relationship significant? How do you know?
##    Yes, p-value is <2e-16
#   2.2 In your own words, what does this mean?
##    That maternal age has a significant effect on the number of DNMs passed from the mother to the child, so older mothers tend to pass on more DNMs.

#2.3 Fit a linear regression model to the data using the lm() function for Maternal.
lm( data = dnm_by_parental_age,
    formula = n_paternal_dnm ~ 1 + Father_age) %>%
  summary()
#   1.1 What is the “size” of this relationship?
##    The size of the relationship is 1.35384
#   1.2 In your own words, what does this mean?
##    The coefficient of 1.35384 indicates that as the father's age increases by 1 year, the number of paternal DNMs is expected to increase by around 1.35384.
#   1.3 Does this match what you observed in your plots in step 2.1?
##    Yes, we can observe a positive trend line in the data

#   2.1 Is this relationship significant? How do you know?
##    Yes, p-value is <2e-16
#   2.2 In your own words, what does this mean?
##    That paternal age has a significant effect on the number of DNMs passed from the father to the child, so older fathers tend to pass on more DNMs.

#2.4 Using your results from step 2.3, predict the number of paternal DNMs for a proband with a father who was 50.5 years old at the proband’s time of birth.
#    Record your answer and your work (i.e. how you got to that answer).
y = Number of Paternal de novo mutations + Paternal Age(Test_Age)
y = 10.32632 + 1.35384(50.5)
y = 78.69524

#2.5 Plot the distribution of maternal DNMs per proband (as a histogram). 
#   In the same panel (i.e. the same set of axes) plot the distribution of paternal DNMs per proband.
#   Make sure to make the histograms semi-transparent so you can see both distributions.
ggplot(data = dnm_by_parental_age) + 
  geom_histogram(aes(x = n_paternal_dnm), fill = "blue", alpha = 0.5) + 
  geom_histogram(aes(x = n_maternal_dnm), fill = "red", alpha = 0.5) + 
  labs(x = "Number of DNMs per proband", y = "Frequency") +
  ggtitle("Distribution of Maternal and Paternal DNMs per Proband")

#2.6 
# test whether there is a significant difference between the number of maternally vs. paternally inherited DNMs per proband.
t_test_result <- t.test(dnm_by_parental_age$n_paternal_dnm, dnm_by_parental_age$n_maternal_dnm, paired = TRUE)
print(t_test_result)

#What would be an appropriate statistical model to test this relationship?
# What statistical test did you choose? Why?
Individually, we can see that the maternal distribution is normal and on the lower end of the x-axis.
While the paternal distribution is also normal but focused on the upper end of the x-axis.
So we can perform a paired t-test on the maternal dnm and the paternal dnm as they separated and have little overlap
# Was your test result statistically significant?
Yes, p-value was less than 0.05 (p-value < 2.2e-16) so there is a significant difference between the number of maternally and paternally inherited DNMs.
This would suggest that the number of DNMs differs between paternal and maternal inheritance, which could imply different factors influencing DNMs from each parent


