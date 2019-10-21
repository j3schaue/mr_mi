# This file was created to attempt to work through the functions
# in imputation_functions.R on a dataset with artificially
# missing data.
# 
# (It's probably going to have a lot of notes to myself and
# not make much sense, haha, ignore me!)

# Source the functions file:
source("sandbox/imputation_functions.R")

# Load libraries (to delete some data):
library(mice) # for multiple imputation 
library(metafor) # for meta-regression

# Extract meta-analytic dataset from metafor:
bd <- metafor::dat.bangertdrowns2004

# Select effect sizes, variances, and some 
# categorical moderator variables:
bd <- bd %>% 
  select(yi, vi, feedback) %>%
  na.omit() # I want to artificially generate missingness, 
# so I'm throwing out data that actually were empirically missing
# (because it's impossible to identify the pattern/reason)

# Make a missingness pattern matrix.
# This has i rows (where i = number of patterns) and j columns
# (where j = number of columns in the dataset/data frame).
# The columns correspond to the columns in the data frame;
# the first two columns are the effect sizes and sampling
# variances, and so on.
# Since I'm just messing with this to understand it,
# I'll just make one pattern with missing data on one predictor
# (let's say 'feedback').
# NOTE: 1s here represent that data is NOT MISSING. 0s here imply
# that it IS missing.

miss_pattern <- matrix(c(1, 1, 0), nrow = 1, byrow = TRUE)
miss_pattern

# Use the ampute() function to delete variables:
bd_mis_obj <- mice::ampute(bd, # data to be amputed
                     prop = 0.5, # proportion of missingness we want
                     patterns = miss_pattern, # missingness pattern
                     mech = "MCAR") # data missing are MCAR

# Extract the data frame with missingness
# and look at it. Note that about 50% of the 'feedback' observations
# are now 'NA', which is correct.
bd_mis <- bd_mis_obj$amp; bd_mis

# Get ambitious and try the impute function!
impute_x_cat(bd_mis)

# Okay nope haha. Too ambitious.

# Ah I think I need to set x equal to the categorical predictor
# with missingness.
x <- bd_mis$feedback
impute_x_cat(bd_mis)

# Oh and y and v to the effect sizes and sampling variances:
y <- bd_mis$yi
v <- bd_mis$vi
impute_x_cat(bd_mis)

# Okay gonna break it down and go step by step (like I probably
# should have done from the beginning, ha):

# Set df (function name for data frame) equal to bd_mis:
df <- bd_mis
# Here I change the labels for df so they match what the code expects:
colnames(df) <- c("y", "v", "x")
# Create "complete"_df by dropping rows with missing data:
comp_df <- df %>% drop_na()
# Create "missing"_df by joining the full dataset (with missing values)
# with the complete dataset and retaining the rows that appear in the
# dataset on the left but not in the dataset on the right
# (which here happens to be the rows missing a covariate).
mis_df <- df %>% setdiff(comp_df)
# ^ Note, not sure why dplyr throws a warning there.

# Create a variable that contains the row numbers for those
# studies with missing data on the specific covariate:
inds <- which(is.na(df$x))
# Working with the COMPLETE data, select every column
# EXCEPT the effect sizes and sampling variances:
comp_df_x <- comp_df %>% 
  select(-y, -v) # That is, creating a matrix of COMPLETE predictors.

# Count the number of predictors (here, 3):
nx <- ncol(comp_df_x)

# This next bit creates dummy-coded variables FROM
# variables that are character or factor types.
# In this specific case, these already are dummy-coded
# variables (each has two possible categories and is represented
# by one variable, which takes on the value of zero or one).
# Therefore, this next bit isn't necessary here (i don't think?)
comp_df_x <- comp_df_x %>% 
  fastDummies::dummy_cols(remove_first_dummy = TRUE)

# Which means that there's no real need to get an updated
# column count here, so I'll just set this equal to 'nx':
ncx <- nx

# And then selecting only dummy-coded predictors should
# just select everything, but to make sure:
comp_df_x <- comp_df_x[, (nx + 1):ncx] %>%
  as.matrix()
# OK that didn't work because it expects nx here to be 4 and
# ncx to be 3, because dummy-coding would drop the first column.
# But no because now that's just [,3:3] which is one column,
# the last column here.
# Actually now I'm confused how it would work anyway. Say there
# actually were 4 binary factors, then it would be [, (4):3]
# wouldn't it?
# Anyhow forget it because this is already a matrix of only dummy-coded
# predictors, hahah.
# Just occurred to me, I think that's because this is only meant
# for one predictor right now, so it makes sense ^^
comp_df_x <- as.matrix(comp_df_x); comp_df_x

post <- get_posterior_stan(t = comp_df$y, 
                          x = comp_df_x, 
                          v = comp_df$v, 
                          fixed = FALSE) # to match the Mr. LL default
# It did the thing! Apparently I have Stan? Haha
# Although it just now occured to me I think this code matches
# one categorical predictor, so it may not have known how to handle
# the others. Going to go back up and drop them.