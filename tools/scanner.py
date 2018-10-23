import numpy as np
import emcee


def make_pred(c,data):
    pred = np.zeros(len(data))
    for ci in range(0, len(coefficients)):
        basis_c_sm_name = coefficients[ci] + "_sm"
        linear_contrib = (c[ci]*self.pred_basis[basis_c_sm_name])
        pred += linear_contrib
    pred = pred + preds['SM']
    return pred
    
def lnprob(c, data, icov):
    pred = make_pred(c, data)
    diff = pred - data
    #print" c = " +str(c) + " pred = " + str(pred) + " data = " + str(data) + " diff = " + str(diff) + " ll = " + str(-np.dot(diff,np.dot(invcov,diff))/2.0)
    return -np.dot(diff,np.dot(self.icov,diff))/2.0

class scanner():
    def init(self,data,icov,preds,coefficients,max_coeff_power):
        #print "initialising scanner"
        self.data = data
        self.icov = icov
        self.preds = preds
        self.coefficients = coefficients
        self.max_coeff_power = max_coeff_power
        pred_basis = {}
        
        #make 'basis predictions'
        for c in coefficients:
            ci_m_name = c + "-"
            ci_p_name = c + "+"
            sigma_sm_ci = np.subtract(preds[ci_p_name], preds[ci_m_name])
            basis_c_sm_name = c + "_sm"
            pred_basis[basis_c_sm_name] = sigma_sm_ci
        self.pred_basis = pred_basis
            
    def make_pred(c):
        pred = np.zeros(len(data))
        for ci in range(0, len(coefficients)):
            basis_c_sm_name = coefficients[ci] + "_sm"
            linear_contrib = (c[ci]*self.pred_basis[basis_c_sm_name])
            pred += linear_contrib
        pred = pred + preds['SM']
        return pred
    
    def lnprob(c, data, icov):
        pred = self.make_pred(c)
        diff = pred - data
        #print" c = " +str(c) + " pred = " + str(pred) + " data = " + str(data) + " diff = " + str(diff) + " ll = " + str(-np.dot(diff,np.dot(invcov,diff))/2.0)
        return -np.dot(diff,np.dot(self.icov,diff))/2.0

    def eval(self, coeff):
        pred =  self.preds['SM']
        chi_sq = 0.0
        for c in range(0, len(self.coefficients)):
            basis_c_sm_name = self.coefficients[c] + "_sm"
            pred += (coeff*self.pred_basis[basis_c_sm_name])
        
        for datum_i in range(0, len(self.data)):
            for datum_j in range(0, len(self.data)):
                chi_sq += (self.data[datum_i] - pred[datum_j])*(self.data[datum_i] - pred[datum_j])*(self.invcov[datum_i][datum_j])
        
        #print "evalutaing point  ------   " + str(self.coefficients[0]) + " = " + str(scenario[0]) + " " + str(self.coefficients[1]) + " = " + str(scenario[1]) + " chi2 = " + str(chi_sq)
        return chi_sq
