

"""
SetRepetitionLabel(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Adds repetition label to raw data object
"""
function SetRepetitionLabel(rawData,params)
    numPar = params[:numPar]
    numSet = params[:numSet]
    
    # Here I want to loop over the profiles and add the "repeition" labels, need to loop over parittions... and adc split 
    i_rep = UInt16(0)
    for i_profile=1:length(rawData.profiles)
        rawData.profiles[i_profile].head.idx.repetition = i_rep
        rawData.profiles[i_profile].head.idx.slice = rawData.profiles[i_profile].head.idx.kspace_encode_step_2
        # Repetitions
        if mod(i_profile,numPar*numSet) == 0
            i_rep += UInt16(1)
        end
    end

    return rawData
end

"""
    FormatRawData(RawData::RawAcquisitionData;single_rep::bool=false,rep_format::Int))

Rearanges the raw data, into [readouts*set,partitions,channels,repetitions]
"""
function FormatRawData(rawData,params;single_rep::Bool=false,rep_format::Int)
    numRead = params[:numRead]
    numInterl = params[:numInterl]
    numPar = params[:numPar]
    numSet = params[:numSet]
    numCha = params[:numCha] 
    numRep = params[:numRep]

    # ToDo: Add functionality for multi-interleaves
    if single_rep == true
        i_prof = Int(((rep_format-1)*numPar*numInterl*numSet)+1)
        numRep = 1
        tmp = Array{ComplexF32}(undef,Int(numRead/numSet),numSet,numInterl,numPar,numCha,1)
        kdata = Array{ComplexF32}(undef,numRead*numInterl,numPar,numCha,1)
    else
        i_prof = 1
        tmp = Array{ComplexF32}(undef,Int(numRead/numSet),numSet,numInterl,numPar,numCha,numRep)
        kdata = Array{ComplexF32}(undef,numRead*numInterl,numPar,numCha,numRep)
    end

    for i_rep=1:numRep
        for i_par=1:numPar
            for i_interl=1:numInterl
                for i_set=1:numSet
                    tmp[:,i_set,i_interl,i_par,:,i_rep] = rawData.profiles[i_prof].data # Original
                    # tmp[:,i_set,i_interl,i_par,:,i_rep] = rawData.profiles[i_prof].data[1:numRead,:]
                    i_prof += 1
                end
            end
        end
    end

    kdata .= reshape(tmp,(:,numPar,numCha,numRep))
    return kdata
end

"""
MergeTrajectory(rawData::RawAcquisitionData,ks_traj::....;skope::Bool=false)

Merges trajectory with 
"""
function MergeTrajectory(rawData,params)
    numSet = params[:set]
    numPar = params[:sl]

    # Here I want to loop over the profiles and add the "repeition" labels, need to loop over parittions... and adc split 
    i_rep = UInt16(0)
    i_par = UInt16(0)
    for i_profile=1:length(rawData.profiles)
        rawData.profiles[i_profile].head.idx.repetition = i_rep
        rawData.profiles[i_profile].head.idx.slice = rawData.profiles[i_profile].head.idx.kspace_encode_step_2
        # Repetitions
        if mod(i_profile,numPar*numSet) == 0
            global i_rep += UInt16(1)
        end
    end

    return rawData
end