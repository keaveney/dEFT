import numpy as np

class scanner():
    def init(self,data,invcov,preds,coefficients,max_coeff_power):
        #print "initialising scanner"
        self.data = data
        self.invcov = invcov
        self.preds = preds
        self.coefficients = coefficients
        self.max_coeff_power = max_coeff_power
        
        pred_basis = {}
    
        #isolate SM--c_i term
        # SM--c_i interference term
        for c in coefficients:
            ci_m_name = c + "-"
            ci_p_name = c + "+"
            sigma_sm_ci = np.subtract(preds[ci_p_name], preds[ci_m_name])
            basis_c_sm_name = c + "_sm"
            pred_basis[basis_c_sm_name] = sigma_sm_ci
        self.pred_basis = pred_basis

    def eval(self, scenario,):
        pred =  self.preds['SM']
        for c in range(0, len(self.coefficients)):
            basis_c_sm_name = self.coefficients[c] + "_sm"
            pred += (scenario[c]*self.pred_basis[basis_c_sm_name])
        
        for datum_i in range(0, len(self.data)):
            for datum_j in range(0, len(self.data)):
                chi_sq = (self.data[datum_i] - pred[datum_j])*(self.data[datum_i] - pred[datum_j])*(self.invcov[datum_i][datum_j])
        
        print "evalutaing point  ------   " + str(self.coefficients[0]) + " = " + str(scenario[0]) + " " + str(self.coefficients[1]) + " = " + str(scenario[1]) + " chi2 = " + str(chi_sq)

    
        return chi_sq
