import numpy as np
class predBuilder:
        def init(self,preds,coefficients,max_coeff_power,c_i_benchmark,cross_terms,inclusive_k_factor):
            self.coefficients = coefficients
            self.predictions = preds
            self.c_i_benchmark = c_i_benchmark
            self.max_coeff_power = max_coeff_power
            self.inclusive_k_factor = inclusive_k_factor
            self.cross_terms = cross_terms
            pred_basis = {}
            #first scale all raw predictions by inclsive k-factor
            unscaled_SM = np.copy(self.predictions['SM'])
            self.predictions['SM'] = np.asarray([x * self.inclusive_k_factor for x in self.predictions['SM']])
            k_term = np.subtract(self.predictions['SM'],unscaled_SM)
            #print "k term  = " + str(self.predictions['SM'])
            
            for c in coefficients:
                ci_m_name = c + "-"
                ci_p_name = c + "+"
                self.predictions[ci_m_name] = np.add(self.predictions[ci_m_name], k_term)
                self.predictions[ci_p_name] = np.add(self.predictions[ci_p_name], k_term)
                #self.predictions[ci_m_name] = np.asarray([x * self.inclusive_k_factor for x in self.predictions[ci_m_name]])
                #self.predictions[ci_p_name] = np.asarray([x * self.inclusive_k_factor for x in self.predictions[ci_p_name]])
            for c in coefficients:
                ci_m_name = c + "-"
                ci_p_name = c + "+"
                # sm--ci interference contribution
                sigma_sm_ci = (np.subtract(self.predictions[ci_p_name], self.predictions[ci_m_name]) / (2*c_i_benchmark))
                # ci--ci squared contribution
                if max_coeff_power > 1.0:
                    #p_add_m = np.add(self.predictions[ci_p_name], self.predictions[ci_m_name])
                    #p_add_m_sub_sm = np.subtract(p_add_m, np.multiply(self.predictions['SM'], 2.0))
                    #sigma_ci_ci = np.divide(p_add_m_sub_sm, np.multiply(self.predictions['SM'],2.0))
                    sigma_ci_ci = np.divide(np.subtract(np.add(self.predictions[ci_p_name], self.predictions[ci_m_name]), np.multiply(self.predictions['SM'], 2.0)), np.multiply(self.predictions['SM'],2.0))
                    #sigma_ci_ci = np.subtract(np.add(self.predictions[ci_p_name], self.predictions[ci_m_name]), (2*self.predictions['SM']) )  / (2*(c_i_benchmark**2))
                    #sigma_ci_ci = (np.add(self.predictions[ci_p_name], self.predictions[ci_m_name]) - (2.0*self.predictions['SM']) ) / (2*(c_i_benchmark**2))
                    # TODO ci--cj interference contribution ("cross-terms")
                    if cross_terms=="true":
                        for d in coefficients:
                            if c != d:
                                ci_cj_m_name = c + "-" + d + "-"
                                ci_cj_p_name = c + "+" + d + "+"
                                sigma_ci_cj=(np.add(self.predictions[ci_cj_m_name],self.predictions[ci_cj_p_name])-(2.0*self.predictions['SM'])-(2.0*(sigma_ci_ci**2))-(2.0*(sigma_cj_cj**2)))/(2*(c_i_benchmark**2))# den. should be 2*(c_i_bench)*(c_j_bench)
                            else:
                                sigma_ci_cj = np.zeros(len(self.predictions['SM']))
                            basis_ci_cj_name = c + "_" + d
                            pred_basis[basis_ci_cj_name] = sigma_ci_cj
                else:
                    sigma_ci_ci = np.zeros(len(self.predictions['SM']))
                basis_ci_sm_name = c + "_sm"
                basis_ci_ci_name = c + "_" + c
                pred_basis[basis_ci_sm_name] = sigma_sm_ci
                pred_basis[basis_ci_ci_name] = sigma_ci_ci
            self.pred_basis = pred_basis
        def make_pred(self,c):
            pred = np.zeros(len(self.predictions['SM']))
            coefficients = list(self.coefficients)
            for ci in range(0, len(coefficients)):
                basis_ci_sm_name = str(coefficients[ci]) + "_sm"
                basis_ci_ci_name = str(coefficients[ci]) + "_" + str(coefficients[ci])
                ci_sm_contrib = (c[ci]*self.pred_basis[basis_ci_sm_name])
                ci_ci_contrib = ((c[ci]**2)*self.pred_basis[basis_ci_ci_name])
                pred += ci_sm_contrib
                if (self.max_coeff_power > 1):
                    pred += ci_ci_contrib
        
                #add option for cross terms (if self.max_coeff_power == ) ...etc.
                #print "ci  in make_pred " + str(self.pred_basis[basis_ci_sm_name])
            pred = pred + self.predictions['SM']
            return pred
