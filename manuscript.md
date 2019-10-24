River discharge can be retrieved from SWOT observations

A data-driven approach led to a plausible estimation of prior
probability distributions

Hydraulic geometry constraints improved river discharge estimation
accuracy across all metrics

Introduction {#sec:orge11834f}
============

Water is essential for all ecosystems and civilizations, and better
understanding of global hydrology benefits numerous fields of study.
Currently, the best available tool to understand the water cycle at any
place on earth is a gauging station, where in situ instruments
(typically) record water level and convert this quantity to mass
flux/discharge via empirically calibrated functions. These gauges are
well understood and accurate when properly maintained, but gauge
maintenance requires physical site access at considerable expense.
Gauges are therefore declining worldwide as they are defunded or as
sharing their information becomes too politically sensitive
[@hannah2011; @gleason2017]. The loss of these primary observations
impairs our ability to tune land surface models, which in turn affects
our understanding of the impacts of climate change, our ability to
manage scarce water resources, and our general confidence in the state
of the global hydrosphere. Since all models are primarily a function of
the quality of their input data, additional primary information is
needed to augment gauges and improve hydrologic understanding worldwide.

Remote sensing is gaining attention as a valuable tool as hydrologists
attempt to fill the temporal and spatial gaps in global observation
networks [@mccabe2017]. Early work built on the legacy of gauges, using
satellite observations of width or water surface elevation in exactly
the same manner as a traditional gauge, and this work has continued to
the present day [e.g. @ashmore2006; @smith2008a; @pavelsky2014].
Moreover, rating curves could also be built from observations of other
quantities from space (as opposed to width/water surface elevation), and
this work too has continued [e.g.
@tarpanelli2013; @dijk2016; @bjerklie2018]. Both of these approaches
rely on some in situ data to develop their empirical equations for
discharge, and are successful in doing so. However, in truly ungauged
regions, training data must come from outside the basin and thus invoke
assumptions of hydrologic transferability that could be problematic
[e.g. @coron2012].

To address this need in ungauged basins and supplement existing in-situ
networks with synoptic measurements, NASA/CNES/CSA/UKSA have been
developing the Surface Water and Ocean Topography (SWOT) mission,
designed to provide simultaneous observations of water surface elevation
(WSE), width, and slope on all global rivers wider than 100m
[@biancamaria2016a]. Many methods designed for estimating river
discharge from SWOT operate under the same basic principle, the
so-called Mass conserved Flow Law Inversion (McFLI). This approach
allows discharge estimation in the absence of any in situ data or
transferability assumptions, and attempts to solve an under-constrained
optimization problem for the discharge needed to produce the hydraulic
conditions observed from remote sensing. A number of algorithms have
been developed that can be considered to follow the mass-conserved flow
inversion approach that range in complexity from Markov Chain Monte
Carlo [@durand2014], Hamiltonian Monte Carlo [@hagemann2017], as well as
variational and Kalman Filter data assimilation schemes
[@andreadis2007; @biancamaria2011; @oubanas2018b]. This approach has
also been successfully deployed using optical [@gleason2018; @feng2019],
but SWOT's synoptic observations should provide the most accurate
discharge estimates with acceptable accuracy for numerous river types
[@pavelsky2014a].

Previous SWOT-based McFLI approaches have shown great promise, but their
implementation and evaluation has shown that there are some issues that
could hinder their accuracy and applicability. Firstly, SWOT-based data
assimilation schemes require a hydraulic model to provide predictions of
discharge that can then be inverted when combined with the SWOT
observations. However, these models have tended to be rather complex and
have either significant data requirements or prodigious computational
expense. Furthermore, even the algorithms that do not require a complex
hydraulic model benefit from prior information (e.g. discharge or river
bathymetry) and when that prior is not accurate or unavailable it
hampers the algorithm's performance [@bonnema2016]. Another issue that
afflicts such estimation algorithms is equifinality which simplistically
refers to the existence of multiple sets of variables (discharge, bed
elevation, roughness etc.) that can lead to the same flow profile as the
one observed making the estimation problem ill-posed [@garambois2015].

One way to ameliorate the equifinality issue is to introduce additional
constraints to the optimal estimation problem. At-a-station hydraulic
geometry (AHG) theory has been confirmed both empirically and
theoretically [@gleason2015] and can be independent from the Manning's
equation on which most of the aforementioned algorithms are based on
[@durand2016]. Essentially, AHG emerged from observing that river flow
width, depth and velocity vary linearly (in log-space) with discharge
for a specific river cross section. @ferguson1986 first confirmed that
these power laws are coincidental to common channel geometric form, and
@dingman2007 proposed a simple mathematical model for channel geometry
that reliably yields AHG. These relationships can be incorporated in the
estimation scheme of a SWOT discharge algorithm in order to reduce the
viable parameter space for the mass-conserved flow inversion.

We developed a data assimilation algorithm, the SWOT Assimilated
DiScharge (SADS), that utilizes a simplified hydraulic scheme to
estimate river discharge from SWOT satellite observations. In contrast
to other data assimilation techniques that have demonstrated the
retrieval of discharge from SWOT-like observations, SADS uses a
hydraulic model with minimal data requirements making it ideal for
application globally, and uses a data-driven approach to derive prior
probability distributions for its parameters. We assessed the
performance of the algorithm within the context of synthetic SWOT
observations by sequentially incorporating hydraulic geometry relations
into the estimation of river discharge.

Algorithm description {#sec:org508f22c}
=====================

The SADS algorithm operates on the set of SWOT observables, i.e. WSE,
width and slope and derives estimates of river discharge and its
associated uncertainty using a data assimilation scheme. The
assimilation scheme involves the "first-guess" estimation of hydraulic
variables by combining a forward model and a set of prior probability
distributions before the assimilation of the SWOT observations. The
objective of the algorithm is the estimation of discharge at each river
reach when SWOT observations become available.

Hydraulic model {#sec:orgcf7575d}
---------------

The forward model in the assimilation scheme is based on the Gradually
Varied Flow (GVF) equations, which describe the steady-state,
non-uniform flow in river channels with gradual variations in water
depth and velocity. The general form of the GVF equation [@chow1955] is
$$\dfrac{dY}{dx} = \dfrac{S_{0} - S_{f}}{1 - Fr^{2}}
\label{eq:gvf}$$ where $Y$ is the average water depth, $x$ is the
longitudinal distance, $S_{0}$ is the channel bed slope, $S_{f}$ is the
water surface slope, and $Fr$ is the Froude number. The latter is given
by $$Fr = \dfrac{Q}{W Y \sqrt{g Y}}
\label{eq:froude}$$ and the water surface slope can be calculated from
the Manning equation $$S_{f} = \dfrac{n^{2} Q^{2}}{W^{2} Y^{10/3}}$$

The integration of the governing equations allows for the calculation of
the longitudinal profiles of the water surface of river flow. Equation
[\[eq:gvf\]](#eq:gvf){reference-type="ref" reference="eq:gvf"} can be
solved using a Runge-Kutta method [@rackauckas2017] for each satellite
overpass of a river, and generate a water surface elevation profile
given characteristics of the river channel geometry, roughness
coefficient, and river discharge. In order to incorporate AHG theory
into the model, channel cross sections are modeled following
@dingman2007 with a form that facilitates the association of SWOT
observables and hydraulic geometry coefficients.

![[\[fig:channel\]]{#fig:channel label="fig:channel"}Simplified
schematic of river channel cross-section geometry (adapted from
@dingman2007).](figures/channel.pdf){width=".9\\linewidth"}

Figure [\[fig:channel\]](#fig:channel){reference-type="ref"
reference="fig:channel"} shows a schematic of the cross section
geometry, where $W^{\ast}$, $W$, $Y_{m}^{\ast}$, $Y_{m}$ are the
bankfull width, water surface width, bankfull depth, and maximum depth
respectively. The latter are related to the average width and flow
depth: $$Y = \left( \frac{r}{r+1} \right) Y_{m}
\label{eq:dingman1}$$
$$W = W^{\ast} \left( \frac{Y_{m}}{Y_{m}^{\ast}} \right)^{1/r}
\label{eq:dingman2}$$ Hydraulic geometry relations can be incorporated
into the SADS GVF model through the relationships between the AHG
coefficients and the channel cross-section geometric variables. These
relationships are derived by rewriting the AHG relations in terms of the
Manning equation and river channel geometry (see @dingman2018 for
details) and are given by
$$W = \alpha Q^{b} = {W^{\ast}}^{5r/3\delta} \left( \frac{Y^{\ast}}{R} \right)^{-5/3\delta} \left(\frac{1}{n} S_{f}^{1/2} \right)^{-1/\delta} Q^{1/\delta}
\label{eq:ahg1}$$
$$Y = c Q^{f} = {W^{\ast}}^{-r/\delta} \left( \frac{Y^{\ast}}{R} \right)^{1/\delta} \left(\frac{1}{n} S_{f}^{1/2} \right)^{-r/\delta} Q^{r/\delta}
\label{eq:ahg2}$$
$$U = k Q^{m} = {W^{\ast}}^{-2r/3\delta} \left( \frac{Y^{\ast}}{R} \right)^{2/3\delta} \left( \frac{1}{n} S_{f}^{1/2} \right)^{(1+r)/\delta} Q^{2r/3\delta}
\label{eq:ahg3}$$ where $R = (1 + r) / r$ and $\delta = 1 + (5 / 3) r$.

It becomes apparent from Equations
[\[eq:dingman1\]](#eq:dingman1){reference-type="ref"
reference="eq:dingman1"} - [\[eq:ahg3\]](#eq:ahg3){reference-type="ref"
reference="eq:ahg3"} that the channel geometry parameter $r$ is very
important. The parameter reflects the river channel shape, with $r=1$
corresponding to a triangular channel and as $r$ increases the channel
banks become steeper and the bottom becomes flatter leading to a
rectangular channel for $r \rightarrow \infty$. There are various
methods to estimate $r$ from channel geometry or other measurements
[@moramarco2019] as well as from the AHG coefficients (if they can be
calculated otherwise), but here we treat $r$ as a stochastic variable
with an associated probability distribution.

Data assimilation scheme {#sec:orgc387af5}
------------------------

The assimilation algorithm employed in the implementation of SADS
presented here is the Local Ensemble Transform Kalman Filter (LETKF)
[@hunt2007]. The LETKF is a variant of the Ensemble Kalman Filter
[@evensen2003], that combines a prior probability distribution of state
variables (e.g. river discharge) with direct or indirect observations
(in this case, water surface elevation and width) to generate an optimal
estimate (i.e. analysis). The prior distribution is represented by the
model error covariance $\bm{P^{b}}$ that is calculated from a
$k\mathrm{-size}$ ensemble of model states (background ensemble),
$\bm{x}_{i}^{b}: i = 1,2,...,k$. The observations $\bm{d^{o}}$ and their
uncertainty can be represented by mapping the state variables to
observation space via $\bm{d}=H(\bm{x})$ and an error covariance matrix
$\bm{R}$ respectively. The analysis (i.e. posterior) mean can be
calculated from
$$\overline{\bm{x}}^{a} = \overline{\bm{x}}^{b} + \bm{X^{b}} \bm{\overline{w}^{a}}
\label{eq:letkf1}$$
$$\bm{\overline{w}^{a}} = \bm{P^{a}} \left(\bm{D}^{b}\right)^{T} \bm{R}^{-1} \left(\bm{d}^{o} - \overline{\bm{d}}^{b}\right)
\label{eq:letkf2}$$
$$\bm{P}^{a} = \left[ (k-1) \bm{I} + \left(\bm{D}^{b}\right)^{T} \bm{R}^{-1} \bm{D}^{b}\right]^{-1}
\label{eq:letkf3}$$ where $\overline{\bm{x}}^{b}$,
$\overline{\bm{d}}^{b}$ are the ensemble mean and $\bm{X}^{b}$,
$\bm{D}^{b}$ are the ensemble deviations from the mean for the model
state and model-predicted observations respectively. The $\bm{W}^{a}$
can be used to calculate the analysis ensemble deviations, which can
then be used to reconstruct the entire analysis ensemble, with
$$\bm{W}^{a} = \left[ (k-1) \bm{P}^{a} \right]^{1/2}
\label{eq:letkf4}$$ $$\bm{X}^{a} = \bm{X}^{b} \bm{W}^{a}
\label{eq:letkf5}$$ The LETKF explicitly supports localization by
applying Equations [\[eq:letkf1\]](#eq:letkf1){reference-type="ref"
reference="eq:letkf1"}-[\[eq:letkf5\]](#eq:letkf5){reference-type="ref"
reference="eq:letkf5"} for local patches of the model domain using a
subset of observations for each. In the case of the SADS algorithm, the
model domain consists of a river network with each river being
partitioned into reaches and each reach being partitioned into cross
sections. Therefore the state vector is formed with river discharge at
each river reach, while the truncated (i.e. local) observation vector
consists of the nearest observations in terms of along-channel flow
distance [@garcia-pintado2015]. The observation vector comprises of the
SWOT observations at the cross section locations, and although the
default variable in this study was water surface elevation it can
optionally include width and slope as well.

The estimation of river discharge from observations such as the ones
that SWOT will provide can be difficult when bed elevation and/or
roughness are unknown due to the equifinality issue [@oubanas2018a]. One
approach that can aid in the solution of such problems is regularization
[e.g. @budd2011], wherein additional constraints are introduced in the
form of penalty terms similar to the observation difference
$d^{o} - \overline{d}^{b}$. In the case of river discharge estimation,
the additional constraints can be derived from the at-a-station
hydraulic geometry relations (Equations
[\[eq:ahg1\]](#eq:ahg1){reference-type="ref" reference="eq:ahg1"} -
[\[eq:ahg3\]](#eq:ahg3){reference-type="ref" reference="eq:ahg3"}). In
particular, it can be shown that assimilating "observations" of the form
$W -
\alpha Q^{b} = 0$ for example is equivalent to a Tikhonov regularization
[@johns2008] prior to assimilating the actual observations.

Derivation of priors {#sec:org8ff9670}
--------------------

Ensemble assimilation methods require the definition of a prior
probability distribution from which to generate the ensemble of
background variables [@evensen2004]. Given that our discharge estimation
algorithm needs to be applicable globally, it has to operate on the
assumption of minimal prior knowledge about river discharge as well as
the inputs to the GVF model. Therefore the algorithm starts with
uninformative priors but uses the observations in a data-driven approach
to estimate the necessary prior probability distributions. The inputs to
the GVF model that are not directly observed include discharge, bed
elevation (as well as bed slope), roughness coefficient, and the channel
shape parameter $r$.

We adapted a rejection sampling approach to derive appropriate prior
distributions for the aforementioned variables. Rejection sampling is a
simple technique used to generate samples from a target distribution $T$
using a proposal distribution $P$. Instead of directly sampling from
$T$, the method generates samples from $P$ and accepts/rejects each of
those samples according to likelihood ratio $\dfrac{t(x)}{L p(x)}$ where
$L$ is a constant ($L > 1$), and $t(x)$, $p(x)$ are the density
functions of $T$ and $P$ respectively [@martino2018]. In our case, the
target distribution is the prior distribution of the unobserved variable
(e.g. bed elevation) and the proposal distribution is an uninformative
prior. Since the density function of the target distribution is unknown,
we use the GVF model as a functional to transform both densities $t(x)$
and $p(x)$ to correspond to density functions of water surface elevation
instead of the target variable. The probability density function of WSE
can be estimated from the observations, thus allowing us to calculate
the likelihood ratio and accept/reject the proposed target-variable
value for its prior distribution.

Figure [\[fig:rejection\]](#fig:rejection){reference-type="ref"
reference="fig:rejection"} shows an example of the adapted rejection
sampling approach for estimating the prior distribution of bed
elevation. The algorithm starts with an uninformative prior as the
proposal distribution, from which a set of bed elevation values are
sampled (Figure [\[fig:rejection\]](#fig:rejection){reference-type="ref"
reference="fig:rejection"}a). The uninformative priors are set as
uniform distributions with the bounds for each unobserved variable
described in Table [\[tab:bounds\]](#tab:bounds){reference-type="ref"
reference="tab:bounds"}. Each of the sampled bed elevation values are
used as inputs to the GVF model to simulate an ensemble of WSE values,
with each value in that ensemble corresponding to a bed elevation value.
Using kernel density estimators for the PDFs of the observed and model
WSE, the likelihood ratio defined above can be calculated and each
ensemble value pair can be accepted or rejected (Figure
[\[fig:rejection\]](#fig:rejection){reference-type="ref"
reference="fig:rejection"}b). Subsequently a new distribution can be
calculated from the accepted samples of bed elevation, forming the prior
to be used in the assimilation.

  **Variable**          **Distribution**   **Parameters**
  --------------------- ------------------ -----------------------------------------------------------------------------------
  Bed elevation         Uniform            \[$H_{min}-20$, $H_{min}$\]
  Discharge             Uniform            \[$Q_{mean}/10$, $Q_{mean}*10$\]
  Roughness             Uniform            \[0.01, 0.07\]
  Shape parameter $r$   Estimated          See section [3](#sec:org7ff2e55){reference-type="ref" reference="sec:org7ff2e55"}

  : [\[tab:bounds\]]{#tab:bounds label="tab:bounds"}Distribution type
  and parameters for the uninformative priors used in the rejection
  sampling approach.

![[\[fig:rejection\]]{#fig:rejection label="fig:rejection"}(a) Proposal
PDF for downstream bed elevation and sampled points (showing both the
subsequently rejected and accepted). (b) Estimated PDF (Model) from
stochastic model simulation of upstream water surface elevation with bed
elevation sampled from proposal PDF. Accepted/rejected water surface
elevation points did/didn't fall under the PDF derived from observations
(Observed), with each point corresponding to a bed elevation
sample.](figures/rejection_sampling.pdf){width=".9\\linewidth"}

Commonly used rejection sampling methods may suffer from low
acceptability rates making them inefficient when used to estimate
posterior distributions [@vrugt2018]. In the case of the SADS algorithm
we use such methods to identify the prior distribution where the
parameter space can be relatively large as long as it contains
behavioral solutions for the particular river. Moreover, the GVF model
is simple enough to make the generation of large populations to sample
from computationally efficient. Following previous studies that
evaluated the performance of SWOT discharge algorithms, we chose to only
use a mean annual flow ($Q_{mean}$) estimated from a global water
balance model as a prior estimation of mean flow [@durand2016] for the
purposes of consistency. In cases where additional information on river
discharge or bed elevation is available, it can be used to theoretically
improve the performance of the SADS algorithm.

Geomorphological classification {#sec:org7ff2e55}
===============================

The at-a-station hydraulic geometry relations that were implicitly
incorporated into the GVF model via the derivations of @dingman2007
depend on the specification of the shape parameter $r$. In the absence
of simultaneous data on discharge and width/depth, the estimation of $r$
becomes difficult and therefore the approach of treating it
stochastically within the SADS assimilation algorithm becomes necessary.
Nonetheless, an approach that would constrain the probable values of $r$
for each river would theoretically improve the accuracy of the discharge
estimation.

We developed a two-pronged approach to constrain the $r$ parameter with
the first part of the approach involving a global predictive model for
$r$ using SWOT-observable variables as predictors, while the second
component pertaining to a supervised expert classification framework
from river planform geometry. Both the estimation and classification
were trained on empirically derived $r$ values from HYDRoSWOT
[@canova2016], a collection of over 200,000 USGS field hydraulic
measurements (originally made for calibrating rating curves for USGS
gages) across the continental US and Alaska. Observed $r$ values equate
to $f/b$, where $f$ and $b$ are the exponents of depth and width
at-a-station hydraulic geometry relationships, respectively, fit at each
cross section. Stations with physically impossible AHG results, such as
$f < 0$ or $f > 1$ were filtered out, as were stations with less than 20
measurements.

The unsupervised predictive model uses a random forest algorithm to
derive predicted $r$ values from river widths, bankfull widths, and
water surface elevations: all observable by SWOT. A large number (500)
of random regression trees were run, with one variable randomly sampled
at each split. Using HYDRoSWOT data, WSE was defined as the sum of mean
depth and datum height at-a-station, and bankfull width was defined as
the maximum width of a given station's inter-quartile range of width
measurements. This model successfully explains 96% of variance in the
dataset's $r$ values and was successfully validated on an independent
set of the data ($R^{2}=0.98$, RMSE=0.23). Note that the model predicts
a range of $r$ values at-a-station, despite our definition of $r$ as
constant at-a-station. This is a function of there being a range of
predictor variables yet a single $r$ at each station, but in the case of
the SADS algorithm it allows us to construct distributions of potential
$r$ mean values for a given cross section.

![[\[fig:rivers\]]{#fig:rivers label="fig:rivers"}Classes of river
planform geometry defined by the expert classification
algorithm.](figures/riv_schematic.pdf){width=".9\\linewidth"}

The supervised expert classification scheme relies on assessing river
planform geometry and making an informed decision of reach geomorphology
(Figure [\[fig:rivers\]](#fig:rivers){reference-type="ref"
reference="fig:rivers"}). First, $r$ was conceptually mapped onto
classic river classification frameworks
[@church2006; @rosgen1994; @schumm1985]. Given the hydraulic definition
of $r$, highly stable and meandering to straight rivers should exhibit
large $r$ values, while unstable and frequently avulsing multi-threaded
systems will exhibit low $r$ values as they move more rapidly in
vertical space than horizontal (as meandering systems do). HYDRoSWOT
empirical $r$ values were qualitatively assessed to determine class
thresholds. It was noted that sinuous single channels generally had $r$
values between 1 and 10, straight "canal"-like channels had values over
10, and $r$ values below 1 resembled both braided channels and
ostensibly wide, shallow channels near dams (verified visually in Google
Earth). It was decided that for the purposes of this algorithm, $r$
values below 1 corresponded to just braided reaches. The distinction
between anastomosing and braiding is inherently subjective, though in
keeping with the nomenclature of @durand2016, these multi-threaded
systems are simply referred to as braided. Class thresholds were finally
set according to these observations and then compared against the
frequency distribution of HYDRoSWOT $r$ values as a qualitative check.
The expert framework's class bounds approximate the following bins:
0-50th percentile, 50-75th percentile, 75th to outlier threshold, and
the outliers. Within SADS the bounds for each river class were used to
truncate the probability distribution of the $r$ parameter, with the
mean and variance of that distribution (assumed to be Gaussian) derived
from the unsupervised predictive model. Table
[\[tab:rivers\]](#tab:rivers){reference-type="ref"
reference="tab:rivers"} show the class of each river from the 18 case
studies using the geomorphological classification framework, along with
the bounds for the $r$ parameter for each class.

  **Classification**                      **Rivers**                                                                                                 
  --------------------------------------- ---------------------------------------------------------------------------------------------------------- --------------
  Somewhat Sinuous Single-Channel Reach   GaronneDownstream, GaronneUpstream, Kanawha, MississippiDownstream, Ohio, Po, SacramentoUpstream, Severn    1 \< r \< 5
  More Sinuous Single-Channel Reach       Cumberland, MississippiUpstream, SacramentoDownstream, Seine, Wabash                                        5 \< r \< 10
  Straight/Canal Reach                    StLawrenceDownstream, StLawrenceUpstream                                                                      r \> 10
  Braided                                 Ganges, Platte, Tanana                                                                                         r \< 1

  : [\[tab:rivers\]]{#tab:rivers label="tab:rivers"}Classification of
  each river from the 18 case studies using the geomorphological
  framework.

Experimental design {#sec:orgaab81e3}
===================

We performed four sets of experiments, each one intended to evaluate the
impact of hydraulic geometry constraints and formulations on the
assimilation of SWOT satellite observations. In each experiment,
different configurations of the SADS algorithm were implemented and
evaluated in terms of river discharge estimation. These experiments
included: 1) the simplified assumption of rectangular river channel
cross sections; 2) using channel cross sections with an empirical but
generic $r$ parameter; 3) deriving the $r$ parameter utilizing our
geomorphological classification scheme; 4) the inclusion of Equations
[\[eq:ahg1\]](#eq:ahg1){reference-type="ref" reference="eq:ahg1"} -
[\[eq:ahg3\]](#eq:ahg3){reference-type="ref" reference="eq:ahg3"} as
observational constraints in the assimilation.

We used the dataset from @durand2016, which includes the output from
hydrodynamic model simulations of 18 rivers, to evaluate the different
configurations of the SADS algorithm. Synthetic SWOT observations were
generated from the output of the hydrodynamic simulations that included
daily values of spatially-variable (i.e. longitudinal profiles) water
surface elevation, slope, and top width corresponding to different
flows. These data were simulated after forcing the different models with
a prescribed channel bathymetry, model-derived or observed inflows and
water elevation as upstream and downstream boundary conditions
respectively. Validation results of these models showed that they
realistically reproduced water surface elevation and inundation extent
[e.g. @jung2012], making them appropriate for the purpose of studies
evaluating river discharge from SWOT-like satellite observations. The
synthetic SWOT observations were derived by adding zero-mean Gaussian
errors (temporally uncorrelated) to the observables with standard
deviations of 5 cm, 5 m, and 0.1 cm/km for water surface elevation,
width, and slope respectively. Prior probability distributions were
estimated for discharge, bed elevation, channel roughness and geometry
parameter ($r$) using the SWOT synthetic observations and methods
described in sections [2.3](#sec:org8ff9670){reference-type="ref"
reference="sec:org8ff9670"} and
[3](#sec:org7ff2e55){reference-type="ref" reference="sec:org7ff2e55"}.
Other inputs to the GVF model included the bankfull width and depth
which were calculated as the maximum observed surface water width and
elevation respectively after subtracting the assumed bed elevation
(derived from the prior PDF that is estimated) for the latter.

The performance of the different configuration of the SADS algorithm
were evaluated with five metrics including relative root mean squared
error (RRMSE), normalized root mean squared error (NRMSE), relative bias
(rBias), Kling-Gupta efficiency (KGE) and Nash-Sutcliffe efficiency
(NSE). RRMSE, given by
$\sqrt{\dfrac{1}{N}\sum\limits_{i=1}^{N}{\left( \dfrac{Q_{i} - \hat{Q_{i}}}{Q_{i}}
\right)^{2}}}$, NRMSE, given by
$\sqrt{\dfrac{1}{N}\sum\limits_{i=1}^{N}{\left(
\dfrac{Q_{i} - \hat{Q_{i}}}{\overline{Q}} \right)^{2}}}$, and rBias,
given by
$\dfrac{1}{N}\sum\limits_{i=1}^{N}{\left( \dfrac{Q_{i} - \hat{Q_{i}}}{Q_{i}}
\right)}$, can help characterize the deviation and mean of the model
residuals and despite their shortcomings [e.g. @ehret2011] have been
some of the most widely used error metrics [@dawson2007]. NSE, given by
$1 -
\sum\limits_{i=1}^{N}{\left(Q_{i} - \hat{Q_{i}}\right)^{2}} \big/
\sum\limits_{i=1}^{N}{\left(Q_{i} - \overline{Q}\right)^{2}}$ is another
widely used goodness-of-fit metric that quantifies the portion of the
observed variance explained by the model. Finally, the KGE presents an
interesting decomposition of the NSE that facilitates the analysis of
different components of the hydrologic signal and is given by
$1-\sqrt{(\rho - 1)^{2} + \left(
\dfrac{\mu_{Q}}{\hat{\mu_{Q}}} - 1 \right)^{2} + \left(
\dfrac{\sigma_{Q}}{\hat{\sigma_{Q}}} - 1\right)^{2}}$ where $\rho$ is
the Pearson correlation between observations and model predictions,
$\mu_{Q}$ and $\sigma_{Q}$ are the mean and standard deviation of the
observations, and $\hat{\mu_{Q}}$ and $\hat{\sigma_{Q}}$ are the mean
and standard deviation of the model predictions [@gupta2009a].

Results {#sec:org4a86246}
=======

Prior distributions {#sec:org1163ca9}
-------------------

Given the importance of prior information on the performance of any
Bayesian schemes [@hagemann2017], it is worthwhile examining whether the
method of deriving the prior PDFs that is part of the SADS algorithm can
reproduce the observed distribution of river discharge. Note that this
prior PDF derivation method is part of all four configurations of the
SADS algorithm implemented and evaluated in our study. Apart from
visually comparing the two distributions, we can quantitatively assess
their similarity using an $f\mathrm{-divergence}$ metric. Here we choose
the Hellinger distance that is bounded between 0 and 1 and can be
considered the probabilistic equivalent of the commonly used Euclidean
distance [@cam2000]. The Hellinger distance between two probability
distributions with densities $f(x)$ and $g(x)$ can be expressed as
$H^{2} = 1 - \int{\sqrt{f(x) g(x)} dx}$. If we assume that the two
probability distributions are lognormal, a good approximation in the
case of river discharge, $H$ can be practically calculated from
$$H^{2} = 1 - \sqrt{\frac{2 \sigma_{1} \sigma_{2}}{\sigma_{1}^{2} + \sigma_{2}^{2}}} exp \left[- \frac{(\mu_{1} - \mu_{2})^{2}}{4 (\sigma_{1}^{2} + \sigma_{2}^{2})} \right]$$
where $\mu$ and $\sigma$ are the means and standard deviations of the
two log-transformed distributions.

![[\[fig:priors\]]{#fig:priors label="fig:priors"}Empirical and
estimated prior distributions of river discharge (normalized by each
river's Truth maximum values for display purposes) for each river case
study. Also shown is the Hellinger distance metric between the two
distributions of each river (smaller values indicate higher degree of
similarity).](figures/sads_priors.pdf){width=".9\\linewidth"}

Figure [\[fig:priors\]](#fig:priors){reference-type="ref"
reference="fig:priors"} shows the comparison of the empirical PDF,
derived from the "true" discharge, and the estimated prior PDF used in
the SADS assimilation algorithm. The comparison is shown in the form of
a set of violin plots that correspond to the probability densities
grouped for each river. The discharge values have been normalized by
each river's "true" maximum value in order to display the PDFs of all
the case studies. With the exception of the Upstream Sacramento River,
SADS appears to reproduce the range of the prior distribution of
discharge. The mode and shape, including dispersion and skewness, were
captured well in the cases of the Ganges, Upstream Garonne, Ohio, Po,
Severn, and Wabash Rivers. On the other hand, the priors for the rest of
the case studies were not able to reproduce some of the properties of
the observed distribution as well potentially affecting the efficacy of
the data assimilation algorithm.

Nonetheless, results from the SADS prior PDF estimation resulted in a
better constraint for the subsequent assimilation. The similarity
between the estimated prior distribution and the corresponding truth is
confirmed in Figure [\[fig:priors\]](#fig:priors){reference-type="ref"
reference="fig:priors"} that also shows the Hellinger distance metric
for each river apart from the direct PDF comparison. The computed
Hellinger distances ranged from 0.002 (Po River) to 0.209 (Platte River)
and a median value of 0.031 suggesting that the SWOT observations and
the rejection sampling method led to an informative prior.

Rectangular channel {#sec:org419405c}
-------------------

Discharge estimation approaches from remote sensing observations that
employ hydraulic models globally utilize the assumption of a rectangular
channel in the absence of other data [e.g. @yoon2012]. Therefore, the
first experiment we performed involved the same approximation for the
channel cross sections, i.e. setting $r$ to a very high value (1e6), in
the GVF model. From the geomorphological classification, most of the
rivers do not appear to be well represented with a rectangular channel.
Nonetheless, the SWOT observations could contain adequate information to
compensate for the uncertainty in the forward model. Figure
[\[fig:rect\]](#fig:rect){reference-type="ref" reference="fig:rect"}
shows hydrographs of the SADS-estimated discharge when using the
rectangular channel approximation compared to the true discharge for
each river case study. In some cases the assimilation is able to
reproduce discharge with reasonable accuracy despite the simplicity in
channel geometry. In particular, the Cumberland and the Seine showed
relatively high efficiencies with NSE of 0.69 and bias of -16%. The St
Lawrence Upstream and Mississippi River reaches also showed relatively
good performance with RMSEs between 24.3 and 41.3%. For the rest of the
cases, the performance of SADS with a rectangular channel ranged from
capturing the temporal variability with some bias (e.g. the Po, Severn
and Sacramento) to completely missing the discharge dynamics (e.g. the
downstream Garonne and St Lawrence Rivers). The largest relative bias in
the rectangular channel estimates are found for the Tanana, Wabash,
Severn, upstream Sacramento, and Ganges with values ranging from -224%
to -604%. Overall, only 6/18 rivers had a positive NSE ranging from 0.12
to 0.69, 12/18 had a positive KGE ranging from 0.06 to 0.76, and 5/18
had NRMSE or RRMSE below 50% ranging from 19% to 48%. The rivers with
the worst negative NSEs were the Wabash (-31.7), the upstream Sacramento
(-26.1), the Tanana (-24.4), and the downstream St Lawrence (-213.8).
Interestingly, the Tanana had a positive KGE value which was reflected
on the RMSEs that were better than the other "worst-NSEs" rivers (73%
versus 157% to 376% in terms of NRMSE). The median values for the error
metrics were 0.11 for the KGE, -0.46 for the NSE, 75% for the NRMSE, 97%
for the RRMSE, and -83% for the relative bias. The relatively poor
results suggest the inappropriateness of the rectangular channel
approximation, although the algorithm might be successfully compensating
for that in some cases by estimating an effective bed elevation and
channel roughness.

![[\[fig:rect\]]{#fig:rect label="fig:rect"}Hydrographs from the SADS
algorithm when assuming a rectangular cross section compared to true
discharge for each of the 18 case study
rivers.](figures/sads_rectangular.pdf){width="\\textwidth"}

AHG channel {#sec:org03f335b}
-----------

In order to improve on the river channel representation in the GVF
forward model, we used the HYDRoSWOT dataset to estimate a more
realistic $r$ shape parameter. The second experiment that we performed
derived a distribution of $r$ from the entire HYDRoSWOT dataset, thus
ignoring any information on the classification of each river. From the
histogram of the data, it was apparent that a log-normal distribution
could be fit and the resulting parameters were $\mu=0.95$ and
$\sigma=1.02$ for the prior probability distribution of $r$. Figure
[\[fig:sads\_r\]](#fig:sads_r){reference-type="ref"
reference="fig:sads_r"} shows time series of the SADS-estimated
discharge when using a river channel with an $r$ parameter sampled from
the aforementioned generic log-normal distribution, along with the true
hydrographs.

The generic $r$ parameter improves the discharge estimation accuracy for
13/18 rivers with the exceptions being the Cumberland, Platte,
downstream Sacramento, Tanana, and upstream St Lawrence. For the cases
where there was improvement, the median KGE increased from 0.06 to 0.33,
the NSE from -0.42 to 0.51, and the NRMSE/RRMSE decreased from 85 and
107% to 43 and 71% respectively. In contrast, the performance metrics
worsened for the five other rivers with KGE decreasing from 0.31 to
-0.06, NSE from -1.02 to -1.54, and RRMSE increasing from 73 to 80% in
terms of median values. Somewhat surprisingly, the NRMSE actually
decreased for the 5 rivers with worse performance (66 to 57%). In the
case of the Cumberland, it appears that the generic $r$ configuration of
SADS does not perform during the first half of the simulation period
when river slope was steeper than about 0.1 cm/km [@durand2016].
Moreover, this configuration of the algorithm does not estimate the flow
peaks for the downstream Sacramento while also losing the some or all of
the dynamics information for the Platte and upstream St Lawrence. On the
other hand, the volume errors appear to be reduced when a generic $r$
parameter is introduced. For example, the Severn NRMSE is reduced from
158% to 92% while the corresponding numbers for the downstream St
Lawrence were 157% and 36% respectively. Overall, this alternative
approximation for the channel geometry led to an overall improvement
with the median relative bias decreasing to -39%, NRMSE/RRMSE decreasing
to 49/77%, and NSE and KGE increasing to 0.31 and 0.33 respectively.

![[\[fig:sads\_r\]]{#fig:sads_r label="fig:sads_r"}Hydrographs from the
SADS algorithm when using hydraulic geometry relation with a generic $r$
parameter compared to true discharge for each of the 18 case study
rivers.](figures/sads_r.pdf){width="\\textwidth"}

Despite the improvement in discharge estimation accuracy, the generic
$r$ is another approximation similar to the rectangular channel
assumption and does not exploit any of the information from the river's
geomorphology or observations to derive an improved prior probability
distribution for $r$. Utilizing the geomorphological classification
framework we can introduce prior knowledge to aid in the inference of
the unobservable channel geometry. Ideally this geomorphological
information will allow the SADS algorithm to derive a PDF (via the
rejection sampling approach) for the $r$ parameter that has the correct
support and is tighter. Figure
[\[fig:sads\_gr\]](#fig:sads_gr){reference-type="ref"
reference="fig:sads_gr"} shows the comparison between true and
SADS-estimated discharge time series, with the latter using a $r$
parameter distribution estimated from the geomorphological
classification framework.

The geomorphological $r$ (Gr hereafter) configuration improves the
accuracy of the estimated discharge for 13/18 rivers while slightly
decreasing it for 3 rivers. The Seine and the Wabash are the two rivers
where performance is worse across all error metrics. The Gr estimation
over-predicts the peak flows for the Seine river and overestimates the
flow volume for the Wabash despite better capturing the temporal
dynamics (compared to the generic $r$). The largest improvements in
terms of efficiency were the upstream Sacramento, the Platte, and the
Tanana with increases in the NSE of 12.9, 2.3 and 26.5 respectively
(corresponding values for KGE were 1.6, 0.6 and 0.7). The Gr
configuration was able to much better reproduce the temporal dynamics in
both St Lawrence reaches, despite not being able to entirely remove the
bias in these cases. In other rivers where there was improvement, the Gr
hydrograph did not exhibit the same high-frequency variability that the
generic $r$ configuration showed with both the upstream and downstream
Mississippi reaches being the most prominent cases of that behavior. The
Ganges is an interesting case, since the rectangular assumption led to
consistent overestimation of discharge whereas both the generic $r$ and
Gr configurations overestimated the low flows but underestimated the
periods of peak flow agreeing with results from other algorithms
[@durand2016]. Overall the Gr configuration improved the error metrics
compared to both the generic $r$ and rectangular channel approximations,
with the median KGE and NSE being 0.57 and 0.66 while the RMSEs and bias
were reduced to 35%/52% (NRMSE/RRMSE) and -19% respectively.

![[\[fig:sads\_gr\]]{#fig:sads_gr label="fig:sads_gr"}Hydrographs from
the SADS algorithm when using hydraulic geometry relation with a $r$
parameter derived from a geomorphological classification compared to
true discharge for each of the 18 case study
rivers.](figures/sads_Gr.pdf){width="\\textwidth"}

AHG constraints {#sec:org65b2924}
---------------

We performed a fourth experiment to evaluate the impact of applying a
regularization step before the assimilation of SWOT observations. The
regularization explicitly incorporated the AHG power laws into the
discharge estimation, which essentially "tightens" and potentially
shifts the prior distribution of discharge that is being used to
generate the ensemble for the standard assimilation of the SWOT
observations. We used an identical configuration in the forward modeling
with the geomorphological $r$ experiment, therefore the only difference
is the pre-conditioning of the prior discharge ensemble through the
regularization with the AHG relationships. Figure
[\[fig:sads\_ahg\]](#fig:sads_ahg){reference-type="ref"
reference="fig:sads_ahg"} shows hydrographs for all rivers of the true
and SADS-estimated discharge with the latter estimate incorporating the
AHG regularization.

The AHG-constrained configuration (hereafter AHGc) improved the accuracy
of the discharge estimates for 13/18 rivers but degraded performance for
5 rivers. The worse-performing rivers included the Tanana, Ganges,
downstream St Lawrence, Po and the downstream Sacramento with the latter
two showing a very small decrease in efficiencies (KGE from 0.58 to 0.57
for the Po, and 0.43 to 0.40 for the downstream Sacramento). The Tanana
shows a relatively small degradation in terms of KGE but a larger one
for the rest of the error metrics (e.g. NRMSE from 31% to 47% and NSE
from -3.54 to -9.38). The accuracy for the downstream St Lawrence
degrades similarly across all metrics except for NSE, which increased
from -7.47 to -4.87 while showing a modest increase in RMSEs. The Ganges
shows a decrease in accuracy when applying the AHG constraints, but a
comparison of any of the SADS configurations (RRMSE of 115% at best)
with the performance of other SWOT-estimation discharge algorithms
(RRMSE of 52% [@bonnema2016]) suggests that either the AHG relations are
limited in braided rivers or that improvements in the SADS algorithm are
needed.

The improvements in discharge estimates in terms of RMSEs (for the
rivers that did improve) were relatively modest with a median decrease
of 8.6% and 7.8% for NRMSE and RRMSE respectively. The largest
improvement in terms of RMSEs was found for the upstream St Lawrence
with NRMSE and RRMSE being reduced to 9.2% and 8.1% respectively. The
Seine and the Severn also show a large improvement with the AHG
constraints attenuating the overestimated peak flows and better
capturing the low flow periods for both rivers. The NSEs for these two
rivers increase from 0.66 to 0.93, and 0.43 to 0.80 respectively. In the
cases of the two Garonne River reaches, the AHGc is able to better
capture the peak flows as well as reduce the errors for the descending
limbs of the hydrographs (particularly evident for the upstream Garonne
during days 176-279). Overall, for the rivers that did improve with the
AHGc the median increase in NSE and KGE was 0.14 and 0.12 while the
NRMSE, RRMSE and rBias exhibited a median decrease of 8.3%, 6.7% and
7.9% respectively.

![[\[fig:sads\_ahg\]]{#fig:sads_ahg label="fig:sads_ahg"}Hydrographs
from the SADS algorithm when using hydraulic geometry relation with a
$r$ parameter derived from a geomorphological classification and
applying the AHG regularization, compared to true discharge for each of
the 18 case study rivers.](figures/sads_ahg.pdf){width="\\textwidth"}

Examining the overall performance of the different configurations it
becomes clear that the incorporation of hydraulic geometry relations,
via the river channel formulation and the AHG equations, resulted in
improved discharge estimates. Figure
[\[fig:metrics\]](#fig:metrics){reference-type="ref"
reference="fig:metrics"} shows boxplots of the discharge error metrics
from all rivers for the different algorithm configurations: the
rectangular channel approximation (Rect), the generic $r$ parameter (r),
incorporating the geomorphological classification (Gr), and including
the regularization with the AHG constraints (AHG). The median NSE
increases with the addition of constraints and information to the
estimation: -0.46 for Rect, 0.31 for r, 0.66 for Gr, and 0.77 for AHG.
The performance in terms of NSE appears more consistent for Rect than
the generic $r$ with a smaller range of efficiencies despite the median
value being lower. This could be attributed to the variance of the prior
$r$ distributions for each approximation. The range of $r$ values for
the Rect distribution are large and lead to practically an identical
channel shape, while the range of $r$ values in the generic distribution
lead to significantly varying channel shapes affecting the forward model
predictions. The AHG configuration has only 3 rivers with negative NSEs,
compared to 12 for Rect, 6 for r and 4 for Gr.

The same pattern was observed in terms of performance for both the RMSE
metrics, with Ahgc having the smallest median errors of 28.6% (NRMSE)
and 47.9% (RRMSE). Moreover the RMSEs decrease as the incorporation of
hydraulic geometry relations is enhanced. The median values for the
NRMSE (RRMSE) were 74.5% (96.8%) for Rect, 49.4% (77.3%) for r, and
34.7% (52.4%) for Gr. Both the Gr and AHG configurations had a smaller
variance of RMSEs compared with the more generic river channel
approximations, as evidenced by the quartile ranges. The distributions
of RMSEs were relatively skewed for all configurations with some rivers
under-performing for the AHGc. Given that the flow depth (difference
between the water surface and bed elevations) is assimilated essentially
twice (one with Equation [\[eq:ahg2\]](#eq:ahg2){reference-type="ref"
reference="eq:ahg2"} and the other with the SWOT observation), any
errors in the SADS-estimated bed elevation could be the cause of volume
errors being worse for some rivers when the additional AHG constraints
are applied.

The AHGc also outperforms the other three configurations in terms of
both relative bias and KGE. The median bias shows that overall all
configurations are negatively biased, although the AHGc has the smallest
value (-7.4%) while Rect has -83.3%, r has -38.9% and Gr has a rBias of
-19.5%. There are only two rivers where SADS-estimated discharge is
consistently positively biased, the Kanawha and the upstream St Lawrence
although the latter has a relatively small bias of 4.9% for the AHGc.
The KGE metric generally shows that all four configurations have
reasonable performance with consistent ranges of metric values between
the generic river channel approximations. As the application of
hydraulic geometry is enhanced, KGE increases with the median values
being 0.12 for Rect, 0.33 for r, 0.57 for Gr, and 0.63 for AHG.

![[\[fig:metrics\]]{#fig:metrics label="fig:metrics"}Boxplot of
discharge error metrics from all case studies for the different
algorithm configurations, including rectangular channel (Rect), generic
$r$ parameter (r), $r$ parameter derived from the geomorphological
classification (Gr), and same as Gr but with the AHG regularization
applied. Range of y-axis values has been cut off for display
purposes.](figures/sads_metrics.pdf){width="\\textwidth"}

Although the geomorphological classification overall improved the
discharge estimated from the assimilation, there could be cases when a
misclassification could limit the effectiveness of the assimilation
and/or the AHG constraints. The initial impact of a misclassification
would manifest in the erroneous truncation of the prior distribution for
the $r$ parameter, and would propagate in the calculations of the GVF
model leading to a potentially poor estimate of the model-observation
difference and model covariance. We tested this hypothesis by noting
that there were two rivers (Wabash and downstream St Lawrence) where the
SADS algorithm was consistently under-performing despite the increasing
application of the AHG constraints. After changing their classification
and applying SADS with the AHG configuration, we found that the best
results were obtained for $r \in (0, 1]$ which suggests a "braided"
river. Figure [\[fig:slwab\]](#fig:slwab){reference-type="ref"
reference="fig:slwab"} shows hydrographs of the true and SADS-estimated
discharge for the downstream St Lawrence and Wabash when using the
geomorphologically classified (SADS-mcr) and new $r$ (SADS-ccr). The
smaller $r$ bounds led to significantly improved discharge for both
rivers, with efficencies being 0.68 (NSE) and 0.69 (KGE) for the Wabash,
and 0.79 (NSE) and 0.91 (KGE) for the upstream St Lawrence. The
"misclassified" NSEs (KGEs) were -69.7 (-5.33) for the Wabash, and -4.87
(0.01) for the downstream St Lawrence.

Despite the simplicity of the GVF model which makes it difficult to
properly handle braided rivers, they could be approximated as very wide,
single channel rivers. However, if the channel shape parameter is not
commensurate with that approximation the assimilation algorithm would be
hindered due to the erroneous predictions from the forward model. On the
other hand, a reason that the $r<1$ produced better results for the
downstream St Lawrence and the Wabash could be the maximum "observed"
width. The Wabash and downstream St Lawrence have the largest maximum
widths from the rivers that were not classified as braided (11,791 and
15,673 respectively). Since the maximum observed width is used as the
bankfull width in the GVF model, the smaller $r$ could be compensating
for the potential overestimation in channel width (even if the
classification was correct).

![[\[fig:slwab\]]{#fig:slwab label="fig:slwab"}Hydrographs of true and
SADS-estimated discharge with AHG regularization and misclasssified
(SADS-mcr) and correctly classified (SADS-ccr) $r$ parameter for the
downstream St Lawrence (a) and Wabash (b)
Rivers.](figures/stlawrence_wabash.pdf){width="\\textwidth"}

It is worthwhile examining whether the regularization with the AHG
constraints led to estimates of discharge that corresponded to empirical
AHG coefficients convergent with their theoretical values (given by
Equations [\[eq:ahg1\]](#eq:ahg1){reference-type="ref"
reference="eq:ahg1"} - [\[eq:ahg3\]](#eq:ahg3){reference-type="ref"
reference="eq:ahg3"}). We estimated the coefficients via linear
regression of the logarithms of the AGHc discharge and the GVF-predicted
flow depth, width and velocity. Since we used a spatially uniform $r$
channel shape parameter, i.e. the exponents in the AHG equations are
identical along each river, we only calculated and compared the
coefficients at the downstream cross sections. Figure
[\[fig:coeffs\]](#fig:coeffs){reference-type="ref"
reference="fig:coeffs"} shows the estimated downstream AHG exponents
($b$, $f$ and $m$) and coefficients ($\alpha$, $c$ and $k$) compared
with their theoretical values. The width AHG relationship appears to
agree best in terms of the theoretical values of the exponent $b$ and
coefficient $\alpha$ with an $R^{2}$ of 0.98 and 0.99 respectively. This
agreement is reasonable given that width is directly observed and
indirectly used in the assimilation via the forward model predictions.
In contrast, the agreement is rather poor for the flow depth exponent
$f$ and coefficient $c$, with a $R^{2}$ of 0.14 and 0.18 respectively.
The bed elevation values that are estimated from SADS are being used
both in the regularization and the assimilation steps, and any errors in
that estimation could lead to divergent results in the $Q-y$ regression
from the theoretical values. Furthermore, the regularization that was
applied essentially amounts to the modification of the discharge prior
ensemble to weakly impose the AHG constraints [@sugiura2013] and
therefore the assimilation could have resulted in a discharge solution
that is not consistent with the preceding applied AHG constraint. The
velocity exponent $m$ and coefficient $k$ appear to have a modest
agreement with their theoretical values, with a $R^{2}$ of 0.51 and 0.56
respectively, which is somewhat expected since velocity depends on both
flow depth and width (via the Manning's equation).

![[\[fig:coeffs\]]{#fig:coeffs label="fig:coeffs"}Estimated hydraulic
geometry exponents ($b$, $f$, $m$ for width, depth, and velocity
respectively) and coefficents ($\alpha$, $c$, $k$ for width, depth, and
velocity respectively) compared with theoretical values (Equations
[\[eq:ahg1\]](#eq:ahg1){reference-type="ref" reference="eq:ahg1"} -
[\[eq:ahg3\]](#eq:ahg3){reference-type="ref" reference="eq:ahg3"}).
Subscript $D$ indicates theoretical
values.](figures/sads_coeffs.pdf){width="\\textwidth"}

Discussion {#sec:orgdcec72e}
==========

While there are a number of algorithms that have been developed to
estimate river discharge from SWOT observations, their implementation
and evaluation has revealed some issues that motivated the approach
presented in this study. The issues that in some cases made discharge
estimation problematic revolved around the requirement of accurate prior
information, the feasibility of applying the algorithm globally, and
most importantly the equifinality of the hydraulic parameters. Although
SADS is build around a data assimilation scheme, it also uses a
data-driven approach to estimate the prior probability distributions
needed for the assimilation of the SWOT observations. In addition, the
at-a-station hydraulic geometry relations that have been confirmed both
theoretically and empirically [@gleason2015] were incorporated into the
algorithm to ameliorate the parameter equifinality issue by further
constraining the estimation problem. The four sets of experiments, i.e.
algorithm configurations, that we performed to assess the value added by
the AHG relationships involved progressively adding each set of
constraints starting from a baseline rectangular channel assumption. The
configuration that explicitly incorporated the AHG power law equations
and the geomorphological class information outperformed the other three
configurations (rectangular channel, empirically-derived but generic
channel, and channel shape based on the geomorphological
classification). Discharge estimation accuracy was quantified with a
suite of metrics including NSE, KGE, RRMSE, NRMSE and rBias. The median
NSEs for the different configurations were -0.46, 0.31, 0.66 and 0.77,
while the NRMSE was 74.5%, 49.4%, 34.7% and 28.6% respectively with
results being similar in terms of the other three metrics. When taking
into account the misclassification of two rivers from the 18 case
studies, the SADS configurations with the full set of AHG constrains had
only positive NSEs with the sole exception being the braided Tanana
river. The applicability of the simple GVF model for braided rivers is
limited and could be the reason why SADS did not perform as well in the
case of the Tanana river, suggesting that some modifications are needed
to improve the forward modeling. The channel formulation with the $r$
shape parameter allows for flexible representation in the forward model
and facilitates the incorporation of hydraulic geometry constraints and
although it can be estimated with minimal prior information the
algorithm's accuracy could be increased if data are available to
calibrate the different model parameters.

We contend that the results above represent excellent performance for a
discharge algorithm that does not invoke any in-situ data. In order to
contextualize the development and implementation of the SADS algorithm,
we can interface with the discussion found in @durand2016, which
compared the performance of five discharge algorithms. Our present
manuscript uses the same data as that study, and we have also used the
same assumptions and prior data in order to enable direct comparison,
i.e. the only input apart from the SWOT observations was a model-derived
mean annual flow used as a prior estimate of mean flow. From those five
algorithms at least one estimated discharge within 35% RRMSE for each of
the rivers common to our two studies. However, no single algorithm
compared by @durand2016 was able to produce results that performed
consistently across all rivers. In contrast, SADS performs consistently
well across a suite of metrics described above, and is also much less
susceptible to outliers than the algorithms in @durand2016.

Our goal here is not to assess whether SADS outperformed other SWOT
discharge estimation algorithms, as @durand2016 was the first ever
effort to compare these algorithms and systematically evaluate their
performance. The lessons learned from that seminal paper have directly
influenced the development of SADS. The configuration of the SADS
algorithm that included the entirety of the hydraulic geometry
constraints led to a positive NSE for all rivers except for the Tanana
River, so it is worth examining potential reasons for the consistent
performance across varying river environments and flow regimes. We
principally observed from the other algorithm results that synergy among
them could prove fruitful. That is, combination of the best elements of
each algorithm could combine elegantly to yield more consistent
performance. Thus, we were inspired by the hydraulic geometry of AMHG
(later to become BAM) and MFG, the Manning's inversion of GaMO and
MetroMan, the data assimilation in GaMO, and the Markov chain of
MetroMan and BAM. The different components of the SADS algorithm as
implemented here have a direct correspondence with these elements from
each of the other SWOT discharge algorithms. @durand2016 suggested that
the addition of a-priori information should improve the effectiveness of
each algorithm, while a hybrid or synergistic approach could lead to a
SWOT discharge product viable for rivers globally. Our data-driven
derivation of the prior probability distributions was able to derive
informative priors from just the SWOT observations and a mean annual
flow, following directly from this suggestion. Consequently, the
conclusions of @durand2016 appear to be supported by the results of the
implementation and evaluation of the SADS algorithm, and we have
achieved results that are both accurate and consistent across rivers.

Conclusions {#sec:org7f9864b}
===========

The upcoming SWOT satellite mission will offer a unique opportunity to
map river discharge at an unprecedented spatial resolution globally from
its observations of water surface elevation, width and slope. Since
river discharge will be indirectly observed from SWOT, a number of
algorithms have been developed with varying complexity to estimate
discharge from SWOT observables. During the implementation and
evaluation of these algorithms some issues arose that motivated the
approach presented in this study (SADS). Data assimilation has a long
history of successfully been applied with SWOT (or SWOT-like)
observations, and therefore was used as the base to build a new
algorithm. SADS was developed by combining a physically-based and
data-driven approach to estimate the prior probability distributions
needed for the assimilation of the SWOT observations. The data-driven
components of SADS were based on a rejection sampling approach as well
as a geomorphological classification framework, while hydraulic geometry
relationships were incorporated into the algorithm to ameliorate issues
of parameter equifinality by further constraining the estimation
problem. The comprehensive dataset first used in @durand2016 allowed for
evaluating the integration of hydraulic geometry relations with data
assimilation of SWOT observations. A set of four experiments that
progressively added the AHG constraints onto a simple rectangular
channel assumption were performed, with the configuration that
incorporated the full set of AHG constraints giving the most consistent
results in terms of discharge estimation as evidenced when examining
five error metrics.

The results from the development of SADS and its evaluation against a
comprehensive albeit synthetic dataset were encouraging but a number of
limitations need to be addressed through future work. Recent work has
utilized airborne measurements from AirSWOT, a SWOT-like airborne
Ka-band radar, to estimate hydraulic variables
[@altenau2017; @tuozzolo2019]. These datasets could potentially be used
to test and evaluate SADS in order to complement the analysis presented
here. Furthermore, an extension to the classic AHG theory is the
At-many-stations hydraulic geometry (AMHG) showed that the empirical
parameters of AHG (valid at cross-sections) are functionally related
along a river [@gleason2015a]. The incorporation of AMHG could lead to
further constraining the prior and posterior distributions of river
discharge as it is reconcilable with AHG [@brinkerhoff2019], making it a
promising approach for an algorithm such as SADS. Finally, the
implementation of the SADS algorithm presented here is not applicable to
an entire river network, since it does not take into account any flows
at junctions or from lateral tributaries. Given that the objective of
this study was the evaluation of hydraulic geometry constraints on the
assimilation of SWOT satellite observations, the aspects of an
operational implementation of the algorithm were considered to be beyond
the scope of this study. Nonetheless, there are a number of approaches
that could be used to enable the application of SADS over river networks
[e.g. @zhu2011].

Funding for this work was provided by the NASA SWOT Science Team and
Terrestrial Hydrology programs. The authors would like to thank Renato
Frasson and Michael Durand for providing the dataset used in this study.
Data and software code used to generate the results presented herein can
be found at
<https://figshare.com/articles/SWOT_Assimilated_DiScharge/10032311>.
