function create_observer(x0, dynamics, outputfun, params)
    nx = 7   # Dimension of state
    nu = 3   # Dimension of input
    ny = 2   # Dimension of measurements
    Ns = 500 # Number of particles
    dg = MvNormal(Diagonal([3e4,0.001]))          # Measurement noise Distribution
    df = MvNormal(Diagonal([0, 0.0000000001, 0.000000000001, 0.00000000001, 0.0000000001, 0.002296, 0.002296]))          # Dynamics noise Distribution
    dx0 = MvNormal(x0,0.00000001^2*I)   # Initial state Distribution
    vecvec_to_mat(x) = copy(reduce(hcat, x)') # Helper function
    #pf = ParticleFilter(Ns, dynamics, outputfun, df, dg, dx0, p=params)
    pf = ParticleFilter(Ns, dynamics, outputfun, df, dg, dx0, p=params)
    return pf
end

function create_observer_GalOx(x0, dynamics, outputfun, params)
    nx = 7   # Dimension of state
    nu = 3   # Dimension of input
    ny = 2   # Dimension of measurements
    Ns = 500 # Number of particles
    dg = MvNormal(Diagonal([3e7,0.00013]))          # Measurement noise Distribution
    df = MvNormal(Diagonal([0, 0.0000000001, 0.000000000001, 0.00000000001, 0.0000000001, 0.006297, 0.006297]))          # Dynamics noise Distribution
    dx0 = MvNormal(x0,0.00000001^2*I)   # Initial state Distribution
    vecvec_to_mat(x) = copy(reduce(hcat, x)') # Helper function
    #pf = ParticleFilter(Ns, dynamics, outputfun, df, dg, dx0, p=params)
    pf = ParticleFilter(Ns, dynamics, outputfun, df, dg, dx0, p=params)
    return pf
end

function FLUMO_discrete(x,u,p,t; tS=0.001)
    V, I, N, A, D, k_offset1, k_offset2  = x
    n, P0, F0, beta1, b_n, a_n, b_a, a_a = p
    c_D = D ./ V
    k_n = maximum([a_n .* (1 .+ abs(c_D)) .^ b_n .+ k_offset1, 0])
    k_a = maximum([a_a .* (1 .+ abs(c_D)) .^ b_a .+ k_offset2, 0])
    dI = -k_n .* I .- k_a .* I .^ n
    dN = k_n .* I
    dA = k_a .* I .^ n
    Vn = x[1] .+ u[3] #Volume with discrete pulse events
    In = x[2] .+ u[2] .+ dI .* tS #Intermediates with discrete pulse events
    Nn = x[3] .+ dN .* tS #Native protein
    An = x[4] .+ dA .* tS #Aggregated protein
    Dn = x[5] .+ u[1] #Denaturant with discrete pulse events
    return [Vn, In, Nn, An, Dn, k_offset1, k_offset2]
end

function FLUMO_discrete_measurement(x,u,p,t)
    V, I, N, _, D, k_offset1,k_offset2  = x
    n, P0, F0, beta1, b_n, a_n, b_a, a_a = p
    c_D = D ./ V
    k_n = maximum([a_n .* (1 .+ abs(c_D)) .^ b_n .+ k_offset1, 0])
    k_a = maximum([a_a .* (1 .+ abs(c_D)) .^ b_a .+ k_offset2, 0])
    dI = -k_n .* (I./V) .- k_a .* (I./V) .^ n
    dAEW = dI ./ beta1
    F = ((I+N)./V)*F0/(P0./V)
    return [F, dAEW]
end

function FLUMO_discrete_GalOx(x,u,p,t; tS=0.001)
    V, I, N, NC, A, D, k_offset1, k_offset2, k_offset3  = x
    n, P0, F0, beta1, b_n, a_n, b_a, a_a = p
    c_D = D ./ V
    k_n = maximum([a_n .* (1 .+ abs(c_D)) .^ b_n .+ k_offset1, 0])
    k_a = maximum([a_a .* (1 .+ abs(c_D)) .^ b_a .+ k_offset2, 0])
    k_nc = maximum([k_offset3, 0])
    dI = -k_n .* I .- k_a .* I .^ n
    dN = k_n .* I .- k_nc .* NC
    dNC = k_nc .* N
    dA = k_a .* I .^ n
    Vn = x[1] .+ u[3] #Volume with discrete pulse events
    In = x[2] .+ u[2] .+ dI .* tS #Intermediates with discrete pulse events
    Nn = x[3] .+ dN .* tS #Native protein
    NCn = x[4] .+ dNC .* tS #Native protein
    An = x[5] .+ dA .* tS #Aggregated protein
    Dn = x[6] .+ u[1] #Denaturant with discrete pulse events
    return [Vn, In, Nn, NCn, An, Dn, k_offset1, k_offset2, k_offset3]
end

function FLUMO_discrete_GalOx_measurement(x,u,p,t)
    V, I, N, NC, A, D, k_offset1,k_offset2, k_offset3  = x
    n, P0, F0, beta1, b_n, a_n, b_a, a_a = p
    c_D = D ./ V
    k_n = maximum([a_n .* (1 .+ abs(c_D)) .^ b_n .+ k_offset1, 0])
    k_a = maximum([a_a .* (1 .+ abs(c_D)) .^ b_a .+ k_offset2, 0])
    dI = -k_n .* (I./V) .- k_a .* (I./V) .^ n
    dAEW = dI ./ beta1
    F = ((I+N)./V)*F0/(P0./V)
    return [F, dAEW]
end