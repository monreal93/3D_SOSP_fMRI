using Pkg
Pkg.activate("/usr/share/5T4/Alejandro/sosp_vaso/sim/")
# Pkg.activate("/usr/share/PhD_2024/sosp_vaso/sim/")

using Revise
using MRIReco, MRISimulation, MAT, NIfTI, ImagePhantoms
using MIRTjim: jim, prompt; jim(:prompt, true)
using Infiltrator
using ImageTransformations
using Statistics
using FLoops
using PaddedViews
using SphericalHarmonicExpansions
# using Plots, ImageView,
using MRIOperators

using FFTW
FFTW.set_provider!("mkl")

include("../recon/functions/fn_addCorrelatedNoise.jl")
include("../recon/functions/fn_calculateGfactor.jl")
include("../recon/functions/fn_save_nii.jl")

phn_sim = 1            # 1=brain, 0=point (psf), for point (PSF), set all to false
cs_sim = true
cs_recon = true
b0_sim = true
b0_recon = true
t2s_sim = false
t2s_recon = false
high_order_recon = false
order_recon =  1            # Recon order (1,2,3)
add_noise = false
gfactor = false
is2d = false
gfactor_replicas = 2
channels = 32            # 0-channels from CS file, >0 less channels (it will be cropped)
# For PSF only.. T2* and b0
psf_t2s = 23e-3     # T2* in s
psf_b0 = 20         # off-resonance in Hz

# Folder and name of sequence to simulate
folder_sim = "simulations_paper"
scan_sim = ["sb_01"]
traj_type = "nom"
# scan_sim = ["sb_01","sb_02","sb_03","sb_04","sb_05","sb_06","sb_07","sb_08","sb_09","sb_10"]

# Folder and name of sensitivity maps and b0 map to use for simulation
folder = "05192023_sv_paper"
scan = "sv_01"
fieldmap = "s01"
path = string("/usr/share/5T4/Alejandro/sosp_vaso/data/",folder)
path_sim = string("/usr/share/5T4/Alejandro/sosp_vaso/data/",folder_sim)

##### Load parameters
params = matread(string(path,"/acq/",scan,"_params.mat"))
params = params["params"]
params["gen"]["n"][3] = params["gen"]["n"][3]*params["gen"]["kz"]
scan_suffix = string(fieldmap,"_",Int(params["gen"]["n"][1]),"_",Int(params["gen"]["n"][2]),"_",Int(params["gen"]["n_ov"][3]*params["gen"]["kz"]))
params["scan_suffix"] = scan_suffix
params["path_sim"] = path_sim
if phn_sim == 0
    gfactor = false
end

# Some constraints
if high_order_recon
    b0_sim = true  
    traj_type = "sk"
end

###### Load volume for simulation
if phn_sim == 1
    phn = matread(string(path,"/acq/fm_",scan,".mat"))
    # phn = matread(string(path,"/acq/fm_","sv_01",".mat"))
    phn = phn["fieldmap"] 
    phn = dropdims(sqrt.(sum(abs.(phn).^2; dims=5));dims=5)
    phn = phn[:,:,:,1]
    phn = convert(Array{ComplexF64,3},phn)
    phn = imresize(phn,(Int(params["gen"]["n"][1]),Int(params["gen"]["n"][2]),Int(params["gen"]["n"][3])))
    # phn_abs = abs.(phn)
    # phn_abs = (phn_abs.- minimum(last,phn_abs))./(maximum(last,phn_abs)-minimum(last,phn_abs))
    # phn_nii = NIVolume(phn_abs)
    # niwrite(string(path_sim,"/sim/reference_phn.nii"),phn_nii)
else
    # phn = zeros(Int(params_pulseq["gen"]["n_ov"][1]),Int(params_pulseq["gen"]["n_ov"][2]),Int(params_pulseq["gen"]["n_ov"][3]))
    # phn[Int(params_pulseq["gen"]["n_ov"][1]/2),Int(params_pulseq["gen"]["n_ov"][2]/2),Int(params_pulseq["gen"]["n_ov"][3]/2)] = 1
    phn = zeros(Int(params["gen"]["n"][1]),Int(params["gen"]["n"][2]),Int(params["gen"]["n"][3]))
    phn[Int(params["gen"]["n"][1]/2),Int(params["gen"]["n"][2]/2),Int(params["gen"]["n"][3]/2)] = 1
    phn = Array{ComplexF64}(phn)
    # phn_nii = NIVolume(abs.(phn))
    # niwrite(string(path_sim,"/sim/reference_phn_point.nii"),phn_nii)
end
if is2d
    sl = Int(params["gen"]["n"][3]/2)
    phn = phn[:,:,sl]
end
phn_orig = deepcopy(phn)

###### Load Coil Sensitivities
cs = matread(string(path,"/acq/cs_",scan_suffix,".mat"))
cs = cs["coil_sens"]
cs = convert(Array{ComplexF64,4},cs)
cs = imresize(cs,(Int(params["gen"]["n"][1]),Int(params["gen"]["n"][2]),Int(params["gen"]["n"][3])))
cs_orig = deepcopy(cs)
# cs = reverse(cs,dims = 2)

##### Load B0 map
if b0_sim
    b0 = matread(string(path,"/acq/b0_",scan_suffix,".mat"))
    b0 = b0["b0"]
    # b0 = convert(Array{ComplexF64,3},b0)
    b0 = imresize(b0,(Int(params["gen"]["n"][1]),Int(params["gen"]["n"][2]),Int(params["gen"]["n"][3])))
    
    if phn_sim == 0
        b0 *= 0
        b0[Int(params["gen"]["n"][1]/2),Int(params["gen"]["n"][2]/2),Int(params["gen"]["n"][3]/2)] = psf_b0*2*pi*im
    end

    b0_orig = deepcopy(b0)
end

##### Load t2s map
if t2s_sim
    if phn_sim == 0
        t2s = zeros(size(phn))
        t2s[Int(params["gen"]["n"][1]/2),Int(params["gen"]["n"][2]/2),Int(params["gen"]["n"][3]/2)] += 1/psf_t2s
    else
        t2s = matread(string(path,"/acq/t2s_",scan_suffix,".mat"))
        t2s = t2s["t2str_map"]
        t2s = convert(Array{ComplexF64,3},t2s)
        t2s = imresize(t2s,(Int(params["gen"]["n"][1]),Int(params["gen"]["n"][2]),Int(params["gen"]["n"][3])))
    end

    t2s_orig = deepcopy(t2s)
end

if gfactor
    ##########################  Fully sampled trajectory
    ###### Load simulation parameters
    full_scan = string(scan_sim[1][1:2],"_full")
    params_pulseq_full = matread(string(path_sim,"/acq/",full_scan,"_params.mat"))
    params_pulseq_full = params_pulseq_full["params"]

    ###### Load trajectory
    traj_full = matread(string(path_sim,"/acq/",full_scan,"_ks_traj_nom.mat"))
    traj_full = traj_full["ks_traj"]

    # Normalizing the traj to -0.5 to 0.5
    traj_full = [traj_full["kx"][:] traj_full["ky"][:] traj_full["kz"][:]]
    traj_full = permutedims(traj_full,[2,1])

    kx_fov_full = maximum([abs.(minimum(traj_full[1,:])) abs.(maximum(traj_full[1,:]))])*2
    ky_fov_full = maximum([abs.(minimum(traj_full[2,:])) abs.(maximum(traj_full[2,:]))])*2
    kz_fov_full = maximum([abs.(minimum(traj_full[3,:])) abs.(maximum(traj_full[3,:]))])*2

    traj_full[1,:] = traj_full[1,:]./kx_fov_full
    traj_full[2,:] = traj_full[2,:]./ky_fov_full
    traj_full[3,:] = traj_full[3,:]./kz_fov_full

    # Creating/loading time vector
    times_full = params_pulseq_full["gen"]["t_vector"][:]		    # Times vector for B0 correction

    # if its 3D, repeat the times vector for each slice
    if !is2d
        times_full = repeat(times_full,Int(params_pulseq_full["gen"]["n"][3]/params_pulseq_full["gen"]["kz"]))
        times_full = vec(times_full)
    end

    ks_full = permutedims(traj_full,(2,1))
    ks_full = Trajectory(traj_full,Int64(params_pulseq_full["spi"]["interl"]),Int64(params_pulseq_full["gen"]["ro_samples"]);times = times_full)

    ###### Simulation
    params_full = Dict{Symbol, Any}()
    params_full[:simulation] = "fast"
    params_full[:trajName] = "Spiral"
    params_full[:AQ] = params["gen"]["acqTR"]
    params_full[:times] = times_full

    # acqData_full = simulation(ks_full, phn, correctionMap=[]; senseMaps=cs, params_full)

    ##########################
end

########################## Under-sampled trajectories
for i = 1:Int(length(scan_sim))
# @floop for i in 1:Int(length(scan_sim))
    
    ###### Load simulation parameters
    params_pulseq = matread(string(path_sim,"/acq/",scan_sim[i],"_params.mat"))
    params_pulseq = params_pulseq["params"]

    # Let's do the matrix size square
    params_pulseq["gen"]["n_ov"][1] = params_pulseq["gen"]["n_ov"][2]
    params_pulseq["gen"]["n"][1] = params_pulseq["gen"]["n"][2]

    params_pulseq["gen"]["n"][3] = params_pulseq["gen"]["n"][3]*params_pulseq["gen"]["kz"]

    ###### Load trajectory
    traj = matread(string(path_sim,"/acq/",scan_sim[i],"_ks_traj_nom.mat"))
    traj = traj["ks_traj"]

    if is2d
        traj["kx"] = traj["kx"][:,sl]
        traj["ky"] = traj["ky"][:,sl]
        traj["kz"] = traj["kz"][:,sl]
    end

    # Normalizing the traj to -0.5 to 0.5
    traj = [traj["kx"][:] traj["ky"][:] traj["kz"][:]]
    traj = permutedims(traj,[2,1])

    kx_fov = maximum([abs.(minimum(traj[1,:])) abs.(maximum(traj[1,:]))])*2
    ky_fov = maximum([abs.(minimum(traj[2,:])) abs.(maximum(traj[2,:]))])*2
    kz_fov = maximum([abs.(minimum(traj[3,:])) abs.(maximum(traj[3,:]))])*2

    traj[1,:] = traj[1,:]./kx_fov 
    traj[2,:] = traj[2,:]./ky_fov 
    traj[3,:] = traj[3,:]./kz_fov

    
    # Normalize k-space trajectory and create Trajectory object
    if traj_type == "sk"
        ks_traj = matread(string(path_sim,"/acq/",scan_sim[i],"_ks_traj_sk.mat")); ks_traj = ks_traj["ks_traj"]
        # ToDo: Sync Skope trajectory to raw data...

        # Swap dimension 3 with 4, Spherical harmonics h2=z, Skope data h2=y.
        ks_traj[4,:] , ks_traj[3,:] = ks_traj[3,:],ks_traj[4,:]

        # ToDo: Do I want to repeat them, one per repetition??
        ks_traj_high = ks_traj
        ks_traj_high = repeat(ks_traj_high, outer=(1,Int(params_pulseq["gen"]["n_ov"][3])))
        # Taking only the specified order 
        ks_traj_high = ks_traj_high[1:Int((order_recon+1)^2),:]

        # # Normalizing
        # for i=axes(ks_traj_high,1)
        #     ks_traj_high[i,:] = ks_traj_high[i,:]./(maximum([abs(minimum(ks_traj_high[i,:])),abs(maximum(ks_traj_high[i,:]))])*2)
        # end

        ks_traj = ks_traj[[2,4,3],:]

        # Repeat trajectory to number of partitions
        ks_traj[3,:] = ks_traj[3,:].-params_pulseq["gen"]["n_ov"][3]/2*(params_pulseq["gen"]["del_k"][3])

        tmp = deepcopy(ks_traj)
        for i=1:Int(params_pulseq["gen"]["n_ov"][3]-1)
            tmp[3,:] = tmp[3,:].+(params_pulseq["gen"]["del_k"][3])
            ks_traj = hcat(ks_traj,tmp)       
        end

        # Normalizing
        for i=axes(ks_traj,1)
            ks_traj[i,:] = ks_traj[i,:]./(maximum([abs(minimum(ks_traj[i,:])),abs(maximum(ks_traj[i,:]))])*2)
        end
        
    elseif traj_type == "nom"
        ks_traj = matread(string(path_sim,"/acq/",scan_sim[i],"_ks_traj_nom.mat")); ks_traj = ks_traj["ks_traj"]
        # Normalizing
        ks_traj["kx"] = ks_traj["kx"]./(maximum([abs(minimum(ks_traj["kx"])),abs(maximum(ks_traj["kx"]))])*2)
        ks_traj["ky"] = ks_traj["ky"]./(maximum([abs(minimum(ks_traj["ky"])),abs(maximum(ks_traj["ky"]))])*2).*(-1)
        ks_traj["kz"] = ks_traj["kz"]./(maximum([abs(minimum(ks_traj["kz"])),abs(maximum(ks_traj["kz"]))])*2)
        ks_traj = hcat(ks_traj["kx"][:],ks_traj["ky"][:],ks_traj["kz"][:])
        ks_traj = permutedims(ks_traj,[2,1])
    end

    ks_traj = convert(AbstractMatrix{Float32},ks_traj)
    times = repeat(params_pulseq["gen"]["t_vector"][:],Int(params_pulseq["gen"]["n_ov"][3]))
    times = convert(Vector{Float32},times)

    ks = Trajectory(ks_traj,1,Int(params_pulseq["gen"]["ro_samples"]*params_pulseq["spi"]["interl"]); 
            times=times,TE=params_pulseq["gen"]["TE"],AQ=params_pulseq["gen"]["ro_time"], 
            numSlices=Int(params_pulseq["gen"]["n_ov"][3]),circular=true)

    if is2d
        traj = traj[1:2,:]
    end

    # # Resizing sensitivities and phantom to simulation matrix size
    # global cs = imresize(cs,(Int(params_pulseq["gen"]["n"][1]),Int(params_pulseq["gen"]["n"][2]),Int(params_pulseq["gen"]["n"][3])))
    # global phn = imresize(phn,(Int(params_pulseq["gen"]["n"][1]),Int(params_pulseq["gen"]["n"][2]),Int(params_pulseq["gen"]["n"][3])))
    # global phn_abs = abs.(phn)
    # global phn_abs = (phn_abs.- minimum(last,phn_abs))./(maximum(last,phn_abs)-minimum(last,phn_abs))

    # Pad zeros to match FOV of pulseq data
    mtx_pulseq = params_pulseq["gen"]["n"][1],params_pulseq["gen"]["n"][2],params_pulseq["gen"]["n_ov"][3]*params_pulseq["gen"]["kz"]
    mtx_pulseq = trunc.(Int,mtx_pulseq)
    mtx = params["gen"]["n"][1],params["gen"]["n"][2],params["gen"]["n"][3]
    mtx = trunc.(Int,mtx)
    # mtx = size(cs)[1:3]

    if cs_sim
        cs = deepcopy(cs_orig)
        if is2d
            cs = cs[:,:,sl,:]
            cs = reshape(cs,(size(cs)...,1))
            cs = permutedims(cs, (1,2,4,3))
        end
        if mtx_pulseq[1:2] > mtx[1:2]
            tmp = zeros((mtx_pulseq...,size(cs,4)))
            tmp1 = Int.((mtx_pulseq.-mtx)./2)
            if tmp1[3] == 0; tmp1 = tmp1.+(0,0,1); end
            if tmp1[2] == 0; tmp1 = tmp1.+(0,1,0); end 
        else
            tmp = zeros((mtx...,size(cs,4)))
            tmp1 = (1,1,1)
        end
        if mtx_pulseq[3] > mtx[3]
            tmp = zeros(size(tmp)[1:2]...,mtx_pulseq[3],size(cs,4))
            tmp1 = tmp1.-(0,0,1)
            tmp1 = tmp1.+(0,0,Int((mtx_pulseq[3].-mtx[3])./2))
        end
        tmp = convert(Array{ComplexF32,4},tmp)
        if is2d
            tmp = tmp[:,:,1,:]
            tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1,:] = cs
            tmp = reshape(tmp,(size(tmp)...,1));
            tmp = permutedims(tmp, (1,2,4,3));
        else
            tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1,tmp1[3]:tmp1[3]+mtx[3]-1,:] = cs
        end

        # Downsampling if requested:
        if mtx_pulseq < mtx
            tmp = imresize(tmp, mtx_pulseq)
        end
        if channels != 0
            tmp = tmp[:,:,:,1:Int(floor(size(tmp,4)/channels)):end]
        end

        global cs = deepcopy(tmp)
    end

    phn = deepcopy(phn_orig)
    if mtx_pulseq[1:2] > mtx[1:2]
        tmp = zeros((mtx_pulseq))
        tmp1 = Int.((mtx_pulseq.-mtx)./2)
        if tmp1[3] == 0; tmp1 = tmp1.+(0,0,1); end
        if tmp1[2] == 0; tmp1 = tmp1.+(0,1,0); end
    else
        # tmp = zeros((mtx...,size(cs,4)))
        tmp = zeros(mtx)
        tmp1 = (1,1,1)
    end
    if mtx_pulseq[3] > mtx[3]
        tmp = zeros(size(tmp)[1:2]...,mtx_pulseq[3])
        tmp1 = tmp1.-(0,0,1)
        tmp1 = tmp1.+(0,0,Int((mtx_pulseq[3].-mtx[3])./2))
    end

    tmp = convert(Array{ComplexF32,3},tmp)
    tmp_orig = deepcopy(tmp)

    if is2d
        tmp = tmp[:,:,1]
        tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1] = phn
    else
        tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1,tmp1[3]:tmp1[3]+mtx[3]-1] = phn
    end
    # Downsampling if requested:
    if mtx_pulseq < mtx
        tmp = imresize(tmp, mtx_pulseq)
    end
    global phn = deepcopy(tmp)
    global phn_abs = abs.(phn)
    global phn_abs = (phn_abs.- minimum(last,phn_abs))./(maximum(last,phn_abs)-minimum(last,phn_abs))

    if gfactor
        ###### Simulation of fully sampled data
        acqData_full = simulation(ks_full, phn, correctionMap=[]; senseMaps=cs, params_full)
        # Scaling down raw data, seems like simulation makes it into 0-1, but real scanner data is ~3e-4
        acqData_full.kdata[1] .= acqData_full.kdata[1].*1e-3 
    end

    ###### Simulation of under-sampled data
    params_sim = Dict{Symbol, Any}()
    params_sim[:simulation] = "fast"
    params_sim[:trajName] = "Spiral"
    params_sim[:AQ] = params["gen"]["acqTR"]
    params_sim[:times] = times

    # Do simulation
    @info ("Simulating k-space data ...")
    if b0_sim && t2s_sim
        b0 = deepcopy(b0_orig)
        tmp = deepcopy(tmp_orig)
        if is2d
            b0 = b0[:,:,sl]
            b0 = reshape(b0,(size(b0)...,1))
            b0 = permutedims(b0, (1,2,4,3))
            tmp = tmp[:,:,1]
            tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1] = b0
        else
            tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1,tmp1[3]:tmp1[3]+mtx[3]-1] = b0
        end
        # Downsampling if requested:
        if mtx_pulseq < mtx
            tmp = imresize(tmp, mtx_pulseq)
        end
        global b0 = deepcopy(tmp)
        t2s = deepcopy(t2s_orig)
        tmp = deepcopy(tmp_orig)
        if is2d
            tmp = tmp[:,:,1]
            tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1] = t2s
        else
            tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1,tmp1[3]:tmp1[3]+mtx[3]-1] = t2s
        end
        # Downsampling if requested:
        if mtx_pulseq < mtx
            tmp = imresize(tmp, mtx_pulseq)
        end
        global t2s = deepcopy(tmp)
        acqData = simulation(ks, phn, t2s+b0; senseMaps=cs, params_sim)
    elseif b0_sim
        b0 = deepcopy(b0_orig)
        tmp = deepcopy(tmp_orig)
        # Downsampling if requested:
        if mtx_pulseq < mtx
            # tmp = imresize(tmp, mtx_pulseq)
            b0 = imresize(b0, mtx_pulseq)
        else
            if is2d
                tmp = tmp[:,:,1]
                tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1] = b0
            else
                tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1,tmp1[3]:tmp1[3]+mtx[3]-1] = b0
            end
            b0 = deepcopy(tmp)
        end
        acqData = simulation(ks, phn, b0; senseMaps=cs, params_sim)
    elseif t2s_sim
        t2s = deepcopy(t2s_orig)
        tmp = deepcopy(tmp_orig)
        @infiltrate
        if is2d
            tmp = tmp[:,:,1]
            tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1] = t2s
        else
            tmp[tmp1[1]:tmp1[1]+mtx[1]-1,tmp1[2]:tmp1[2]+mtx[2]-1,tmp1[3]:tmp1[3]+mtx[3]-1] = t2s
        end
        # Downsampling if requested:
        if mtx_pulseq < mtx
            tmp = imresize(tmp, mtx_pulseq)
        end
        global t2s = deepcopy(tmp)
        acqData = simulation(ks, phn, t2s; senseMaps=cs, params_sim)
    elseif cs_sim
        acqData = simulation(ks, phn; senseMaps=cs, params_sim)
    else
        acqData = simulation(ks, phn; senseMaps=[], params_sim)
    end

    # Adding extra parameters to acqData
    acqData.fov =  Tuple(params_pulseq["gen"]["fov"])

    # Scaling down raw data, seems like simulation makes it into 0-1, but real scanner data is ~3e-4
    if phn_sim == 1
        acqData.kdata[1] .= acqData.kdata[1].*1e-3 
    end

    if params_pulseq["gen"]["ro_type"] == "s"
        noise_snr = Float64(5*(params_pulseq["spi"]["interl"]))
        if params_pulseq["spi"]["interl"] > 1
            noise_snr = noise_snr*2
        end
    elseif params_pulseq["gen"]["ro_type"] == "c"
        noise_snr = Float64(2)
    end

    ###### Reconstruction
    params_reco = Dict{Symbol, Any}()

    if is2d
        params_reco[:reconSize] = Int(params_pulseq["gen"]["n"][1]),Int(params_pulseq["gen"]["n"][2])
    else
        # params_reco[:reconSize] = Int(params_pulseq["gen"]["n"][1]),Int(params_pulseq["gen"]["n"][2]),Int(params_pulseq["gen"]["n"][3])
        params_reco[:reconSize] = size(phn)
    end
    params_reco[:regularization] = "L1"
    params_reco[:λ] = 1.e-2
    params_reco[:iterations] = 20
    params_reco[:solver] = "admm"
    params_reco[:method] = "nfft"
    if scan_sim[i][1] == 's'
        params_reco[:rxyz] = params_pulseq["gen"]["kz"]*params_pulseq["spi"]["rxy"]
    elseif scan_sim[i][1] == 'c'
        params_reco[:rxyz] = params_pulseq["gen"]["kz"]*params_pulseq["epi"]["ry"]/params_pulseq["epi"]["pf"]
    end

    # AMM: Temp: Rescaling the B0 map before reconstruction, to have a missmatch btw sim and recon
    @info("Temp.... Reescaling B0 map before recon...")
    @infiltrate
    b0 .= b0*(0.5)

    if cs_recon
        # # Amm: Temp: rotate and reverse CS maps to match simulation data
        # for j=1:size(cs,4)
        #     for i=1:size(cs,3)
        #         cs[:,:,i,j] = rotr90(cs[:,:,i,j]) 
        #     end
        # end
        # cs = reverse(cs; dims=2)
        params_reco[:reco] = "multiCoil"
        params_reco[:senseMaps] = cs
    end
    if high_order_recon
        params_reco[:iterations] = 4
        params_reco[:reco] = "expandedSignalModel"
        params_reco[:higherOrderTrajectory] = convert(Matrix{Float32},ks_traj_high)
        params_reco[:senseMaps] = cs
        params_reco[:correctionMap] = b0
    end
    if b0_recon || (b0_recon && t2s_recon)
        params_reco[:correctionMap] = b0
    elseif t2s_recon
        params_reco[:correctionMap] = t2s
    end

    # load covariance Matrix
    if  params["gen"]["field_strength"] == 7
        cov = matread(string("/usr/share/5T3/Alejandro/sosp_vaso/data/tmp/noise_cov_7T.mat")); cov = cov["C"]
    elseif params["gen"]["field_strength"] == 9
        cov = matread(string("/usr/share/5T3/Alejandro/sosp_vaso/data/tmp/noise_cov_9T.mat")); cov = cov["C"]
    end

    # G-factor map
    if gfactor
        @info ("Calculating G-factor ...")
        params_reco[:path_sim] = params["path_sim"] 
        params_reco[:scan_suffix] = string(scan_sim[i],"_",mtx_pulseq[1],"_",mtx_pulseq[2],"_",mtx_pulseq[3])
        gfactor_map = calculateGfactor(acqData,acqData_full,gfactor_replicas,cov,params_reco)
        gfactor_map = NIVolume(gfactor_map)
        path_gfactor = string(path_sim,"/sim/",scan_sim[i],"_gmap.nii")
        niwrite(path_gfactor,gfactor_map)
    end
    
    ##### Add Correlated noise
    if add_noise && cs_sim
        # @info ("Adding noise to k-space data ...")
        acqData.kdata[1] = addCorrelatedNoise(acqData.kdata[1],noise_snr,cov,Float64(1e1))
    end

    # # #### Temp: trying some direct reconstruction by matrix inversion
    # harmonics = MRIOperators.CalculateSphericalHarmonics((92,94,4),acqData.fov; l =1)
    # E = exp.(im.*3*pi.*transpose(ks_traj_high[2:end,:])*transpose(harmonics[:,2:end]))
    # y = acqData.kdata[1]
    # x = y \ E
    # x = reshape(transpose(x),(92,94,4,4))
    # xx = conj(cs).*x
    # jim(abs.(x))
    # # #####

    @info("stop before recon...")
    @infiltrate

    # Do reconstruction
    @info ("Reconstruction ...")
    Ireco = reconstruction(acqData, params_reco)

    ###### Quality measurments
    # Normalize datasets
    if cs_recon
        Ireco_abs = abs.(Ireco[:,:,:,1,1])
    else
        Ireco_abs = abs.(Ireco[:,:,:,1,:,1])
    end
    Ireco_abs = (Ireco_abs.- minimum(last,Ireco_abs))./(maximum(last,Ireco_abs)-minimum(last,Ireco_abs))
    recon_error = phn_abs.-Ireco_abs

    suffix = "_recon"
    if phn_sim == 0
        suffix = string(suffix,"_point")
    end
    if cs_sim
        suffix = string(suffix,"_cssim")
    end
    if cs_recon
        suffix = string(suffix,"_csrecon")
    end
    if b0_sim
        suffix = string(suffix,"_b0sim")
        if phn_sim == 0
            suffix = string(suffix,Int(psf_b0))
        end
    end
    if b0_recon
        suffix = string(suffix,"_b0recon")
    end
    if t2s_sim
        suffix = string(suffix,"_t2ssim")
        if phn_sim == 0
            suffix = string(suffix,Int(psf_t2s*1e3))
        end
    end
    if t2s_recon
        suffix = string(suffix,"_t2srecon")
    end

    if cs_sim == true && cs_recon ==false
        Ireco = Ireco[:,:,:,1,:,1]
        Ireco = sqrt.(sum((Ireco.^2),dims=4))
        Ireco = Ireco[:,:,:,1]
    end

    # save phn
    phn_abs = abs.(phn)
    phn_abs = (phn_abs.- minimum(last,phn_abs))./(maximum(last,phn_abs)-minimum(last,phn_abs))
    phn_nii = NIVolume(phn_abs; voxel_size=Tuple(params["gen"]["res"].*1e3),time_step=params["gen"]["volTR"])
    if phn_sim == 1
        path_save_recon = string(path_sim,"/sim/",scan_sim[i],"_phn",".nii")
    else phn_sim == 0
        path_save_recon = string(path_sim,"/sim/",scan_sim[i],"_phn_point",".nii")
    end
    niwrite(path_save_recon,phn_nii)

    # save recon
    if phn_sim == 0
        Ireco_nii = NIVolume(real.(Ireco); voxel_size=Tuple(params["gen"]["res"].*1e3),time_step=params["gen"]["volTR"])
    else
        Ireco_nii = NIVolume(abs.(Ireco); voxel_size=Tuple(params["gen"]["res"].*1e3),time_step=params["gen"]["volTR"])
    end
    path_save_recon = string(path_sim,"/sim/",scan_sim[i],"_sim",suffix,".nii")
    niwrite(path_save_recon,Ireco_nii)
    
    # save error
    recon_error = NIVolume(recon_error; voxel_size=Tuple(params["gen"]["res"].*1e3),time_step=params["gen"]["volTR"]);
    path_save_error = string(path_sim,"/sim/",scan_sim[i],"_error",suffix,".nii")
    niwrite(path_save_error,recon_error)
end