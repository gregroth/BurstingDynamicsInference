#calculation of the time scale in interburst survival probability

function timescale_off(model, parameters)
    P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off, PNN, sspTr_off, Pabs = StoThyLiveCell.mo_basics(model, parameters, 40, 1, 1) ;
   
    w = weightsTr_off'*eigvecs(PNN);
    s = inv(eigvecs(PNN))*ones(size(PNN,1));

    ptscales = w'.*s
    tscaletemp = eigvals(PNN)
    tscaletemp[tscaletemp.<=0] .= 1e-10
    tscales = log.(tscaletemp)
    return tscales, ptscales
end


function timescale_off_wosinglet(model, parameters)
    nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PNN, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos = StoThyLiveCell.mo_basics_wosinglet(model, parameters, 40, 1, 1)
   
    w = weightsTr_on_wos'*eigvecs(PNN);
    s = inv(eigvecs(PNN))*ones(size(PNN,1));

    ptscales = real.(w'.*s)
    tscaletemp = real.(eigvals(PNN))
    tscaletemp[tscaletemp.<=0] .= 1e-10
    tscales = log.(tscaletemp)
    return tscales, ptscales
end

function timescale_on(model, parameters)
    P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off, PNN, sspTr_off, Pabs = StoThyLiveCell.mo_basics(model, parameters, 40, 1, 1) ;
   
    weightsAbs_on = ssp[stateAbs_on]./sum(ssp[stateAbs_on])     
    weightsTr_on = weightsAbs_on' * P[stateAbs_on,stateTr_on]./sum(weightsAbs_on' * P[stateAbs_on,stateTr_on])

    PBB = P[stateTr_on,stateTr_on]
    tempdist = weightsTr_on


    w = weightsTr_on*eigvecs(PBB);
    s = inv(eigvecs(PBB))*ones(size(PBB,1));

    ptscales = w'.*s
    tscaletemp = eigvals(PBB)
    tscaletemp[tscaletemp.<=0] .= 1e-10
    tscales = log.(tscaletemp)
    return tscales, ptscales
end


function timescale_on_wosinglet(model, parameters)
    nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PNN, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos = StoThyLiveCell.mo_basics_wosinglet(model, parameters, 40, 1, 1)
   
    PBB = P[stateTr_on,stateTr_on]


    w = weightsTr_on'*eigvecs(PBB);
    s = inv(eigvecs(PBB))*ones(size(PBB,1));

    ptscales = [w[i]*s[i] for i in eachindex(w)]
    tscaletemp = eigvals(PBB)
    ptscales = real.([w[i]*s[i] for i in eachindex(w)]./sum([ptscales[i]*tscaletemp[i] for i in  eachindex(ptscales)]))
    tscaletemp = real.(tscaletemp)
    tscaletemp[tscaletemp.<=0] .= 1e-10
    tscales = log.(tscaletemp)
    return tscales, ptscales
end




function top_timescales(model_th, parameters)
    th_tscaleoff, th_ptscaleoff= timescale_off_wosinglet(model_th, parameters)
        #remove non singl eigen values
        idx = []
        for k in eachindex(th_tscaleoff)
            idx_same = findall(th_tscaleoff.==th_tscaleoff[k])
            if length(idx_same)>1
                prob_tot = th_ptscaleoff[k]
                for j in eachindex(idx_same)
                    if idx_same[j]>k
                        push!(idx,idx_same[j])
                        prob_tot += th_ptscaleoff[idx_same[j]]
                    end
                end
                th_ptscaleoff[k] = prob_tot
            end
        end
        idx = unique(idx)
        deleteat!(th_tscaleoff, idx)
        deleteat!(th_ptscaleoff, idx)
        p =sortperm(th_ptscaleoff)
        th_tscaleoff = th_tscaleoff[p]
        th_ptscaleoff = th_ptscaleoff[p]
        idx_cut = findlast(th_ptscaleoff.<.01)

        toptscale_off = th_tscaleoff[idx_cut+1:end]
        topptscale_off = th_ptscaleoff[idx_cut+1:end]
        top_ts_off = [(topptscale_off[i],abs(toptscale_off[i])) for i in eachindex(topptscale_off)]

    th_tscaleon, th_ptscaleon= timescale_on_wosinglet(model_th, parameters)
        #remove non singl eigen values
        idx = []
        for k in eachindex(th_tscaleon)
            idx_same = findall(th_tscaleon.==th_tscaleon[k])
            if length(idx_same)>1
                prob_tot = th_ptscaleon[k]
                for j in eachindex(idx_same)
                    if idx_same[j]>k
                        push!(idx,idx_same[j])
                        prob_tot += th_ptscaleon[idx_same[j]]
                    end
                end
                th_ptscaleon[k] = prob_tot
            end
        end
        idx = unique(idx)
        deleteat!(th_tscaleon, idx)
        deleteat!(th_ptscaleon, idx)
        p =sortperm(th_ptscaleon)
        th_tscaleon = th_tscaleon[p]
        th_ptscaleon = th_ptscaleon[p]
        idx_cut = findlast(th_ptscaleon.<.01)

        toptscale_on = th_tscaleon[idx_cut+1:end]
        topptscale_on = th_ptscaleon[idx_cut+1:end]
        top_ts_on = [(topptscale_on[i],abs(toptscale_on[i])) for i in eachindex(topptscale_on)]

return top_ts_on, top_ts_off
end