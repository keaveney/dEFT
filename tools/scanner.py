class scanner():
    def init(self,data,invcov,pred_basis):
        print "initialising scanner"
        self.data = data
        self.invcov = invcov
        self.pred_basis = pred_basis
        
    def eval(self, scenario):
        print "evaluating objective function" + str(scenario)
        
        pred =  (self.pred_basis['SM']) + (scenario*self.pred_basis['ctg-'])
    
        for datum_i in range(0, len(self.data)):
            for datum_j in range(0, len(self.data)):
                chi_sq = (self.data[datum_i] - pred[datum_j])*(self.data[datum_i] - pred[datum_j])*(self.invcov[datum_i][datum_j])
    
        return chi_sq
