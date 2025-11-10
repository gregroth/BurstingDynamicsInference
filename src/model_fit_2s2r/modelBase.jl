#set up the topology of the model m_2s2r
#----------------------------------------------------------------------
# Important: provide the name and order of the parameter rates
#
#parameters = (kr1r2,kr2r1,kon1,kon2,koff,kini,delta)
#states = R1ON, R1OFF, R2ON, R2OFF
#----------------------------------------------------------------------



function modelInstance()
    #define the link between the rate matrix and the parameter indices
    Qstate = [0    5    1    0;
            3    0    0    1;
            2    0    0    5;
            0    2    4    0]

    #define the states from where the promoter can initiate
    stateTrlocal = [1,3]
    #define the indices of the initiation rates in the parameter vector; must be the same size than stateTrlocal
    kiniToRate_idx = fill(6,2) #unique kini for the three active states
    #define the number of unique initiation rate
    nbinirate = 1
    #define the index of the degradation rate in the parameter vector
    deltaToRate_idx = 7
    #find the matrix indices where there is a paramter value
    paramToRate_idx = findall(Qstate .>0)
    #define which parameter corresponds to which position in Qstate
    paramToRate_val = Qstate[findall(Qstate .>0)]
    return StoThyLiveCell.StandardStoModel(size(Qstate,1),maximum(Qstate),nbinirate,paramToRate_idx,paramToRate_val,stateTrlocal, kiniToRate_idx, deltaToRate_idx)
end
