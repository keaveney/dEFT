import json
import numpy as np
from numpy.linalg import inv
from tools import yoda2array

class configReader:
    def init(self, filename):
        self.filename = filename
        coefficients = []
        predictions = {}
        with open(filename, 'r') as f:
            config = json.load(f)
        self.params = config
        self.run_name = self.params["config"]["run_name"]
        self.observable = self.params["config"]["data"]["observable"]
        self.bins = self.params["config"]["data"]["bins"]
        self.data = self.params["config"]["data"]["central_values"]
        self.prior_limits = self.params["config"]["model"]["prior_limits"]
        self.tex_labels = self.params["config"]["model"]["tex_labels"]
        self.coefficients = list(self.params["config"]["model"]["prior_limits"].keys())
        self.n_walkers = self.params["config"]["fit"]["n_walkers"]
        self.n_burnin = self.params["config"]["fit"]["n_burnin"]
        self.n_total = self.params["config"]["fit"]["n_total"]
        self.cov = self.params["config"]["data"]["covariance_matrix"]
        self.icov = inv(self.params["config"]["data"]["covariance_matrix"])
        self.x_vals = self.params["config"]["data"]["bins"]
        self.samples = np.asarray(self.params["config"]["model"]["samples"])
        self.kfac = np.asarray(self.params["config"]["model"]["inclusive_k_factor"])

        slabels=[]
        for c in self.coefficients:
            slab = "$" + c + "$"
            slabels.append(slab)
        self.labels = slabels
        x_vals = np.zeros(len(self.x_vals)-1)
        for x_val in range(0, len(self.params["config"]["data"]["bins"])-1):
            x_vals[x_val] = (self.params["config"]["data"]["bins"][x_val] + self.params["config"]["data"]["bins"][x_val+1])/2.0
        self.x_vals = x_vals

        try:
            basestring
        except NameError:
            basestring = str
            
        if(self.params["config"]["model"]["input"] == "numpy"):
            self.predictions = np.asarray(self.params["config"]["model"]["predictions"])
            
        
        if (len(self.predictions) != len(self.samples)):
            raise TypeError('number of samples must equal number of predictions')







        
