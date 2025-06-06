% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/excess.R
\name{ENE}
\alias{ENE}
\title{Estimate Condition-Specific Excess Number of Events in Treatment Group}
\usage{
ENE(
  file,
  reg.model = c("poisson", "NB", "quasipoisson"),
  wt.type = "treat",
  save = TRUE,
  save.plot = TRUE
)
}
\arguments{
\item{file}{path to an Excel file containing the data. The file must contain multiple sheets with the following sheet naming convention: YYYY_Control containing data for control population in year YYYY and YYYY_Treat containing data for population under treatment. Inside each sheet, the condition name/ID must be specified in the first column and the first row contains variable names. Each row contains data for a condition across different agegroups (from column 2 onwards) with the last row contains the population size to be used as offset in the regression model.}

\item{reg.model}{Regression model used to estimate the statistic. The default is Poisson regression with Negative Binomial regression as alternative.}

\item{wt.type}{Which population to use for calculating the 'Overall' Excess estimate as weighted average of the the age-specific Excess estimate. The default is to weigh using the treated (exposed) population. This gives a realistic assessment of the disease burden and spread for the exposed population, and also allows comparisons between the exposed and general populations, accounting for differences in age, so can be interpreted as the effect of exposures if the general population had the same age structure as the population with exposure.}

\item{save}{whether to automatically save the output as a CSV file.}

\item{save.plot}{whether to produce volcano plot and save it as a plotly object in an HTML file.}
}
\value{
A table containing the ENE statistic (per 100,000 person-years) for each condition and agegroup.
}
\description{
This function uses Poisson or Negative Binomial regression to estimate the Excess Number of Events (ENE) per 100,000 person-years. Agegroup-specific and Overall ENE statistic is calculated for each condition.
}
