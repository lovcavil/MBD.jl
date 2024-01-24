using LinearAlgebra
include("../mathfunction_II8.jl")
#include("AppData_door.jl")
export AppData_II8,process_vector,AppDataStruct

struct AppDataStruct
    name::String
    nb::Int
    ngc::Int
    nh::Int
    nc::Int
    NTSDA::Int
    SJDT::Array{Any}
    SMDT::Array
    STSDAT
    q0::Array
    qd0::Array
end

function println0(apps::AppDataStruct)
    println("App Data")
    println("========")
    println("Name: ", apps.name)
    println("Number of Bodies (nb): ", apps.nb)
    println("Number of Geometric Constraints (ngc): ", apps.ngc)
    println("Number of Holonomic Constraints (nh): ", apps.nh)
    println("Number of Constraints (nc): ", apps.nc)
    println("\nSJDT :")
    for row in eachrow(apps.SJDT')
        # println(row)
    end
    println("\nSMDT :")
    println(apps.SMDT)
    println("\nInitial (q0): ", apps.q0)
    println("Initial Derivative (qd0): ", apps.qd0)
    println("========")
end

function process_vector(q::Vector, flags::Vector{Bool}, target::Vector{Vector{Float64}}, func::Function)
    # Ensure `q` and `flags` have the same length
    # if length(q) != length(flags)
    #     error("Vectors 'q' and 'flags' must be of the same length.")
    # end

    # Process elements of `q` based on `flags`
    for (i, flag) in enumerate(flags)
        if flag && i <= length(target)
            func(target[i], q[i])
        end
    end
end

# Function to append `element` of `q` to the `target_vector`
function append_to_vector(target_vector::Vector{Float64}, element)
    push!(target_vector, element)
end

function AppData_II8(app)
    ux = [1, 0, 0]
    uy = [0, 1, 0]
    uz = [0, 0, 1]
    z3 = zeros(3)
    z4 = zeros(4)

    if app == 4   # Rotating Disk with Translating Body

        nb = 2      # Number of bodies
        ngc = 7 * nb  # number of generalized coordinates
        nh = 2      # Number of holonomic constraints
        nhc = 10    # Number of holonomic constraint equations
        nc = nhc + nb  # Number of constraint equations
        nv = ngc - nc
        nu = nc
        NTSDA = 1   # Number of TSDA force elements
    
        # SJDT: Joint Data Table 
        # Format: [t; i; j; sipr; sjpr; d; vxipr; vzipr; vxjpr; vzjpr; a; b or R; mus; mud; mc; nm]
        SJDT = Array{Any}(undef, 28, nh)
        SJDT[:, 1] = Any[4, 1, 0, z3..., z3..., 0, ux..., uz..., ux..., uz..., 0.1, 0.1, 0.08, 0.07, 1, 5]  # Rev - 1 to Ground
        SJDT[:, 2] = Any[5, 1, 2, ux..., z3..., 0, ux..., uy..., ux..., uy..., 0.1, 0.1, 0.08, 0.07, 6, 5]  # Tran. - 1 to 2
    

    
        # SMDT: Mass Data Table (Centroidal with diagonal inertia matrix)
        SMDT = hcat(vcat(10, 10, 10, 10), vcat(5, 5, 5, 5))
    
        # STSDAT: TSDA Data Table
        STSDAT = Array{Any}(undef, 12, NTSDA)
        STSDAT[:, 1] = Any[1, 2, (ux + 10 * uy)..., z3..., 10, 0, 10.1, 0]
    
        # Initial generalized coordinates
        r10 = [0, 0, 0]
        p10 = [1, 0, 0, 0]
        r20 = [1, 0, 0]
        p20 = [1, 0, 0, 0]
        q0 = vcat(r10, p10, r20, p20)
        r1d0 = [0, 0, 0]
        p1d0 = [0, 0, 0, 0]
        r2d0 = [0, 0, 0]
        p2d0 = [0, 0, 0, 0]
        qd0 = vcat(r1d0, p1d0, r2d0, p2d0)
    end
    


    if app == 6  # One Body Translation Along x axis

        nb = 1      # Number of bodies
        ngc = 7 * nb  # number of generalized coordinates
        NTSDA = 0   # Number of TSDA force elements
        nh = 1    # Number of holonomic constraints
        nhc = 5   # Number of holonomic constraint equations
        nc = nhc + nb   # Number of constraint equations
        nv = ngc - nc
        nu = nc   
        SJDT = Array{Any}(undef, 28, nh)
        # SJT(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr;a;b or R;mus;mud;mc;nm];
        # k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
        # 6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
        # si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr;
        # a,b or R=joint dimensions; mus,mud=coefficients of static & dynamic friction.
        # ms=address of first Lag. Mult,; nm= No. of Constr. Eq.

        SJDT[:, 1] = Any[5, 1, 0, z3..., z3..., 0, uz..., ux..., uz..., ux..., 0.1, 0.1, 0.3, 0.25, 1, 4]  # Tran. Bod1-Grnd
    

    
        # SMDT(4,nb): Mass Data Table (Centroidal with diagonal inertia matrix)  
        # SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
        SMDT = [75, 30, 30, 30]
    
        # STSDAT(12,1): TSDA Data Table
        if NTSDA == 0
            STSDAT = zeros(12, NTSDA)
        end
        # STSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F];  
        # T=TSDA No.; i&j=bodies conn.;si&jpr=vectors to Pi&j; K=spring constant;
        # C=damping coefficient; el0=spring free length; F=const. force
    
        # Initial generalized coordinates
        r10 = [0, 0, 0]
        p10 = [1, 0, 0, 0]
        q0 = [r10; p10]
        r1d0 = [1, 0, 0]
        p1d0 = z4
        qd0 = [r1d0; p1d0]
    end


    return nb,ngc,nh,nc,nv,nu,NTSDA,SJDT,SMDT,STSDAT,q0,qd0
end