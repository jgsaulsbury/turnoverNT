% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_spindles.R
\name{plot_spindles}
\alias{plot_spindles}
\title{Plots a set of abundance trajectories as spindles}
\usage{
plot_spindles(occs, ages, plot.ss = TRUE, linesevery = NA, ...)
}
\arguments{
\item{occs}{matrix of the number of observations in each species at each time.
One column for each species, one row for each time slice. Time goes from oldest
at the bottom to youngest at the top.}

\item{ages}{vector containing the ages of each time slice, in years, from
oldest to youngest.}

\item{plot.ss}{boolean indicating whether to plot the per-timeslice sample size
on the right side of the plot.}

\item{linesevery}{value indicating the interval at which to draw horizontal gridlines.
Default is NA (don't draw them).}
}
\value{
plots spindle diagrams depicting the relative abundance of each species
through the timeseries.
}
\description{
A convenient plotting function that takes community composition timeseries in
the same format as that used by fitJ and that produced by simNT. Plots a spindle
for each species, sorted by maximum abundance in the timeseries.
}
\examples{
J <- 10000
tslength <- 500
every <- 5
nsp <- 4
ages <- seq(0,tslength,every)
timeseries <- simNT(startingabs=rep(J/nsp,nsp),ts=ages,ss=1000)
plot_spindles(occs=timeseries$simulation,ages=timeseries$times,linesevery=100)
}
