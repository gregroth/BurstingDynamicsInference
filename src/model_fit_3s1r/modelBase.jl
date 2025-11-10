#set up the topology of the model m_3s1r
#----------------------------------------------------------------------
# Important: provide the name and order of the parameter rates
#
#parameters = (koff, kon, kdown, kup,kini,delta)
#stats: ON, OFF, DeepOFF
#----------------------------------------------------------------------

function modelInstance()
    #define the link between the rate matrix and the parameter indices
    Qstate = [0    1    0;
            2    0    3  ;
            0    4    0 ]

    #define the states from where the promoter can initiate
    stateTrlocal = [1]
    #define the indices of the initiation rates in the parameter vector; must be the same size than stateTrlocal
    kiniToRate_idx = [5] #unique kini for the three active states
    #define the number of unique initiation rate
    nbinirate = 1
    #define the index of the degradation rate in the parameter vector
    deltaToRate_idx = 6
    #find the matrix indices where there is a paramter value
    paramToRate_idx = findall(Qstate .>0)
    #define which parameter corresponds to which position in Qstate
    paramToRate_val = Qstate[findall(Qstate .>0)]
    return StoThyLiveCell.StandardStoModel(size(Qstate,1),maximum(Qstate),nbinirate,paramToRate_idx,paramToRate_val,stateTrlocal, kiniToRate_idx, deltaToRate_idx)
end
