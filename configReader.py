import json
import numpy as np

class configReader:
    def init(self, filename):
        self.filename = filename
        coefficients = []
        predictions = {}

        print "initialising with file " + self.filename
        #Read JSON data into the datastore variable
        with open(filename, 'r') as f:
            config = json.load(f)
        self.params = config
        self.prior_limits = self.params["config"]["model"]["prior_limits"]
        self.coefficients = self.params["config"]["model"]["prior_limits"].keys()
        for predname in self.params["config"]["model"]["predictions"].keys():
            predictions[predname] = np.asarray(self.params["config"]["model"]["predictions"][predname])
        self.predictions = predictions








        
