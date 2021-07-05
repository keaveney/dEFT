# dEFT  - differential Effective Field Theory tool

dEFT is a tool for performing fits of EFT coefficients to HEP data in seconds. 

dEFT employs the Metropolis-Hastings algorithm to efficiently approximate
the likelihood function of the data in the potentially higher-dimensional
space of a typical EFT model. Even with O(10) dimensions, a dEFT fit typically
takes ~30 seconds to run on the desktop.
 
The dEFT philosophy -

dEFT aims to facilitate fast and easy EFT fits on the desktop. A given analysis is entirely defined
by a single json file. A global analysis utlising multiple independent results is performed
by running dEFT with a group of such json files as input. Hence dEFT can be used both to perform
simple, transparent analyses of small datasets and limited sets of EFT operators and complex
analyses of mutiple results each containing many data points and large sets of relevant operators.

Installation  - 

dEFT is run as a python application and has been tested with python 3.X.

dEFT requires a number of packages to be installed that can be easily
obtained with the following pip commands:
 
pip install matplotlib

- Used for plotting. 

pip install numpy

- used for intermediate storage of data and predictions and numerical manipulations
  crucial to the fits. 

pip install emcee
    - implements the Metrolpolis-Hastings method to estimate the N-dimensional likelihood
      function and hence derive the confidence/credible intervals on the EFT coefficients.
      More information on this package can be found here: http://dfm.io/emcee/current/.

pip install corner 
    - used to generate the array of 1- and 2-d scatter plots that visualise the confidence/credible intervals
      More information on this package can be found at https://corner.readthedocs.io/en/latest/install.html.

pip install tqdm
    - used to display progress bars while sampling is running

Building a dEFT analysis:

A dEFT analysis has three basic inputs:

1) The data.
  - This can be as simple as a single data point or as complex
as a multi-differential distribution. If unfolded, differential data
is fitted, the full covariance matrix is highly desirable.

2) The model.
  - The user must specify the EFT coefficients that are to be fitted along with associated
    theoretical predictions for the data. 
  - The predictions must be in the form of a set of "basis", predictions 
 from which predictions for any set of values for the coefficients can be generated. 
 These predictions can be manually written into the json file or a paths  existing yoda
 files containing the predictions may be provided.

3) The fit parameters.
  - The user must specify some extra parameters to control the fit proceedure.
 Some examples are the limits of the fit for each coefficient, the maximum number of
evaluation of the likelihood function and the highest-order terms to be included in
the predictions.   

These three inputs are encapsulated in a single JSON file which entirely
defines a given dEFT analysis.

Once the JSON file has been defined, dEFT can be run with:

python run_dEFT.py analyses/myAnalysis.json 

James Keaveney
