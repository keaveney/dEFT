import numpy as np
class predBuilder:
        def init(self,preds,coefficients,max_coeff_power,c_i_benchmark,cross_terms):
            self.coefficients = coefficients
            self.predictions = preds
            self.c_i_benchmark = c_i_benchmark
            self.max_coeff_power = max_coeff_power

            self.cross_terms = cross_terms
            pred_basis = {}
            for c in coefficients:
                ci_m_name = c + "-"
                ci_p_name = c + "+"
                # sm--ci interference contribution
                sigma_sm_ci = (np.subtract(preds[ci_p_name], preds[ci_m_name]) / (2*c_i_benchmark))
                # ci--ci squared contribution
                if max_coeff_power > 1.0:
                    sigma_ci_ci = (np.add(preds[ci_p_name], preds[ci_m_name]) - (2.0*preds['SM']) ) / (2*(c_i_benchmark**2))
                    # TODO ci--cj interference contribution ("cross-terms")
                    if cross_terms=="true":
                        for d in coefficients:
                            if c != d:
                                ci_cj_m_name = c + "-" + d + "-"
                                ci_cj_p_name = c + "+" + d + "+"
                                sigma_ci_cj=(np.add(preds[ci_cj_m_name],preds[ci_cj_p_name])-(2.0*preds['SM'])-(2.0*(sigma_ci_ci**2))-(2.0*(sigma_cj_cj**2)))/(2*(c_i_benchmark**2))# den. should be 2*(c_i_bench)*(c_j_bench)
                            else:
                                sigma_ci_cj = np.zeros(len(preds['SM']))
                            basis_ci_cj_name = c + "_" + d
                            pred_basis[basis_ci_cj_name] = sigma_ci_cj
                else:
                    sigma_ci_ci = np.zeros(len(preds['SM']))
            
                basis_ci_sm_name = c + "_sm"
                basis_ci_ci_name = c + "_" + c
                pred_basis[basis_ci_sm_name] = sigma_sm_ci
                pred_basis[basis_ci_ci_name] = sigma_ci_ci
            self.pred_basis = pred_basis

        def make_pred(self,c):
            pred = np.zeros(len(self.predictions['SM']))
            
            for ci in range(0, len(self.coefficients)):
                basis_ci_sm_name = self.coefficients[ci] + "_sm"
                basis_ci_ci_name = self.coefficients[ci] + "_" + self.coefficients[ci]
                ci_sm_contrib = (c[ci]*self.pred_basis[basis_ci_sm_name])
                ci_ci_contrib = ((c[ci]**2)*self.pred_basis[basis_ci_ci_name])
                #print "ci_sm_contrib = " + str(ci_sm_contrib)
                pred += ci_sm_contrib
                if (self.max_coeff_power > 1):
                    pred += ci_ci_contrib
                #add option for cross terms (if self.max_coeff_power == ) ...etc.
            pred = pred + self.predictions['SM']
            return pred
