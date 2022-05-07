import pandas as pd
import numpy as np
import math
from scipy.optimize import minimize, rosen, rosen_der

class QuClu:
    def __init__(self) -> None:
        pass


    
    def fun_cu(self,theta, data, k, cl, qq): #cl must be an np.array
        VV = 0
        p = data.shape[1] 
        for i in range(k):
            if np.sum(cl==i)>0:
                nn = np.sum(cl==i)
                xx = data.values[cl==i]
                a  = np.sum((theta+((1-2*theta)*(xx < np.tile(qq[i,], (nn,1))))) * #compute the distance of the entrie dataset wrt each cluster
            np.absolute(xx - np.tile(qq[i,], (nn,1))),axis=1)
                VV =+ np.sum(a)
        n = len(cl)
        VV = VV - p * n * np.log(theta*(1-theta))
        return VV


    def alg_CU(self,data, k = 2,  eps = 1e-8, it_max = 100, B = 30):
        nobs = data.shape[0]
        p = data.shape[1]
        qq = np.zeros((k,p))
        QQ = np.zeros((nobs,k))
        VV_temp = None
        theta=0.5
        VV_list = np.array([math.inf])
        for j in range(B):
            theta = np.random.uniform(0,1)
            for i in range(k):
                qq[i,] = np.apply_along_axis(np.quantile,0,data,((i+1) - 1)/(k - 1) * 0.5 + theta/2) #calculate the quantile of the variable wrt to each clsuter

                #now we calculate the sum(of the different p variable) distances of each points
                QQ[:,i] = np.sum((theta+((1-2*theta)*(data.values < np.tile(qq[i,], (nobs,1))))) *
                   np.absolute(data.values-np.tile(qq[i,], (nobs,1))),axis=1) - p*np.log(theta*(1-theta))
            cl = np.apply_along_axis(np.argmin,1,QQ)
            VV_temp = np.sum(QQ[range(nobs), cl])
            if VV_temp < VV_list[0,] : #check that the new likelihood distance is  less
              VV = VV_temp       #than the one before and  update partition and parameter theta
              cl_true = cl
              theta_true = theta
        # Clustering step
        cl = cl_true
        theta = theta_true
        ratio = 5
        h = 0
        while ratio > eps and h < it_max :
            h = h + 1
            #calculate cluster size
            nk = dict(zip(range(k),np.repeat(0, k, axis=0)))
            values, counts = np.unique(cl, axis=0, return_counts=True)
            nk_new = dict(zip(values,counts))
            nk.update(nk_new)
            for i in range(k):
                if nk[i] > 0 :
                    qq[i,] = np.apply_along_axis(np.quantile,0,data.values[cl==i],theta) #find new vluster quantiles
                    QQ[:,i] = np.sum((theta+((1-2*theta)*(data.values < np.tile(qq[i,], (nobs,1))))) * #compute the distance of the entrie dataset wrt each cluster
                   np.absolute(data.values-np.tile(qq[i,], (nobs,1))),axis=1) - p*np.log(theta*(1-theta))
            #Now update theta

            theta = minimize(self.fun_cu, x0=theta,args = (data, k, cl, qq), bounds=((1e-4,0.99),)).x[0]
            cl = np.apply_along_axis(np.argmin,1,QQ)
            VV_list = np.append(VV_list, np.sum(QQ[range(nobs), cl]))
            ratio = (VV_list[h-1,] - VV_list[h,]) / VV_list[h-1,]
            if h < 5:
                ratio = 2 * eps
        qq1 = pd.DataFrame(qq)
        qq1.columns = data.columns
        res = lt ={ 
            'Vseq' : VV_list,
            'VV' : VV_list[h,],
            'cl' : cl,
            'qq' : qq1,
            'theta' : theta} 
        return res

