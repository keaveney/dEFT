# dEFT  - A desktop Effective Field Theory tool

dEFT is a tool for performing fits of EFT coefficients to HEP data in seconds. 

dEFT employs the Metropolis-Hastings algorithm to efficiently approximate
the likelihood function of the data in the potentially higher-dimensional
space of a typical EFT model. Even with O(10) dimensions, a dEFT fit typically
takes ~30 seconds to run on the desktop.
 
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
