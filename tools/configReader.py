import json
import numpy as np
from numpy.linalg import inv
from tools import yoda2array

class configReader:
    def init(self, filename):
        self.filename = filename
        coefficients = []
        predictions = {}
        #print "initialising with file " + self.filename
        #Read JSON data into the datastore variable
        with open(filename, 'r') as f:
            config = json.load(f)
        self.params = config
        self.run_name = self.params["config"]["run_name"]
        self.observable = self.params["config"]["data"]["observable"]
        self.prior_limits = self.params["config"]["model"]["prior_limits"]
        self.coefficients = self.params["config"]["model"]["prior_limits"].keys()
        self.n_walkers = self.params["config"]["fit"]["n_walkers"]
        self.n_burnin = self.params["config"]["fit"]["n_burnin"]
        self.n_total = self.params["config"]["fit"]["n_total"]
        self.cov = self.params["config"]["data"]["covariance_matrix"]
        self.icov = inv(self.params["config"]["data"]["covariance_matrix"])
        self.x_vals = self.params["config"]["data"]["bins"]
        self.cross_terms = self.params["config"]["model"]["cross_terms"]
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

        for predname in self.params["config"]["model"]["predictions"].keys():
            if (isinstance(self.params["config"]["model"]["predictions"][predname], basestring)):
                predictions[predname] = yoda2array.convert(self.params["config"]["model"]["predictions"][predname], self.params["config"]["model"]["histname"])
            else:
                predictions[predname] = np.asarray(self.params["config"]["model"]["predictions"][predname])
        self.predictions = predictions







        
