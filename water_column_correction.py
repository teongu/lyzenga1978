"""
This module performs the water column correction as described in Lyzenga 1978. 

The implementation is strongly inspired from the work of @jkibele, some functions are directly his work (https://github.com/jkibele/OpticalRS/blob/master/OpticalRS/Lyzenga1978.py) but some modifications have been made. 

You can find an example of how to use this code in the notebook "Lyzenga 1978 demo.ipynb".

In the following documents, `N` stands for the number of bands.

-----
Reference: 
Lyzenga, D.R., 1978. Passive remote sensing techniques for mapping water depth and bottom features. Appl. Opt. 17, 379–383. doi:10.1364/AO.17.000379
"""


import numpy as np
from sklearn.linear_model import LinearRegression

def compute_Xi(bands, list_deep_water):
    """
    Compute the transformed band Xi, as defined in Equation 7 of Lyzenga 1978: Xi=ln(Li-Lsi)
    -----
    Input: 
        - bands : list of 2D array. Each 2D array corresponds to a band of the image that needs to be corrected.
        - list_deep_water : list of float. Each float represents the mean of a given band on a deep water zone. 
    bands and list_deep_water must have the same length.
    -----
    Output:
        - list of 1D array. Each array correspond to the Xi, which has been reshaped in 1D (mandatory to perform the correction)
    """
    if len(bands)!=len(list_deep_water):
        raise Exception("Lists must have the same length.")
    list_Xi=[]
    for band,dw in zip(bands,list_deep_water):
        list_Xi.append(np.log(band.reshape(-1)-dw))
    return(np.array(list_Xi))

def linear_regression(list_X,Z):
    """
    Perform the linear regression as described in Equation B1 of Lyzenga 1978: Xi=ai-bi*z
    -----
    Input:
        - list_X : list of 1D array. Each 1D array corresponds to the points of Xi where the regression has to be performed.
        - Z : 1D array (or list). Points where the depth is known and where the regression has to be performed.
    Z and each Xi in list_X must be the same length.
    -----
    Output:
        - list of the slopes for each band (i.e. [b0,b1,...,bN])
        - list of the intercepts for each band (i.e. [a0,a1,...,aN])
    """
    list_slope=[]
    list_intercept=[]
    for X in list_X:
        if len(X)!=len(Z):
            raise Exception("Xi and Z must be the same length.")
        model = LinearRegression()
        model.fit(Z.reshape(-1,1), X)
        list_slope.append(-model.coef_[0])
        list_intercept.append(model.intercept_)
    return(np.array(list_slope),np.array(list_intercept))


def Aij(b):
    """
    Return the coordinate system rotation parameters described in Equations B5, B6 and B7, 
    that are needed to compute the depth invariant index of Equation 8.
    ----------
    Input:
        - b : 1D array (or list). Slopes from regressing depths against logged radiance values.
    -------
    Output:
        - 2D array. Matrix Aij as defined in Equations B5, B6 and B7 of Lyzenga 1978. 
    """
    N = len(b)
    A = np.empty((N,N), dtype='float')
    for i in range(N):
        for j in range(N):
            if i==N-1:
                # Equation B2
                A[i,j] = b[j]*(b**2).sum()**(-0.5)
            elif j<=i:
                # Equation B5
                A[i,j] = b[i+1]*b[j]*((b[:i+1]**2).sum()**(-0.5))*((b[:i+2]**2).sum()**(-0.5))
            elif j==i+1:
                # Equation B6
                A[i,j] = -1.*((b[:i+1]**2).sum()**(0.5))*((b[:i+2]**2).sum()**(-0.5))
            else:
                A[i,j] = 0.
    return(A)


def Y_i(i,A,X):
    """
    Calculate band i of the coordinate rotation described in Equation 8 of Lyzenga 1978.
    -----
    Input:
        - i : int. The band number of Y to calculate.
        - A : 2D array. The coordinate system rotation parameters Aij as described in Appendix B of Lyzenga 1978.
        - X : 2D array of size (N,len(Xi)). The array of all the 1D Xi. 
    -----
    Output:
        - 1D array, of size (len(Xi)). Band Y_i as described in equation 8 of Lyzenga 1978.
    """
    return((A[i,:] * X.transpose()[...,:]).sum(1))


def depth_invariant(A,X):
    """
    Calculate the N-1 depth-invariant bands and single depth dependent band by coordinate rotation described in equation 8 of Lyzenga 1978.
    Parameters
    -----
    Input:
        - A : 2D array. The coordinate system rotation parameters Aij as described in Appendix B of Lyzenga 1978.
        - X : 2D array of size (N, len(Xi)). The array of all the 1D Xi.
    -----
    Output:
        - 2D array of size (N, len(Xi)). The N-1 depth-invariant bands and single depth dependent band as described in equation 8 of Lyzenga 1978.
    """
    Yi = []
    for i in range(X.shape[0]):
        Yi.append(Y_i(i,A,X))
    return(np.array(Yi))

def reshape_di(di, bands):
    """
    Reshape the depth invariant bands Yi, originally of size (len(Xi),), to the shape of the input images.
    -----
    Input:
        - di : 2D array of size (N, len(Xi)). Depth invariant bands as described in equation 8 of Lyzenga 1978.
        - bands : list of 2D array. Each 2D array corresponds to a band of the image that needs to be corrected.
    Output:
        - list of 2D array. Each 2D array correspond to a band of the depth invariant image.
    """
    reshaped=[]
    for i in range(di.shape[0]):
        reshaped.append(np.reshape(di[i,:],np.shape(bands[i])))
    return(reshaped)