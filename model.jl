using Parameters, Distributions, StatsBase, Random, Match, DataFrames

# health status
@enum HEALTH SUS LAT PRE ASYMP INF REC DED UNDEF

# define an agent and all agent properties
Base.@kwdef mutable struct Student
    idx::Int64 = 0                              # id 
    health::HEALTH = SUS                        # health status
    swap::HEALTH = UNDEF                
    sickfrom::HEALTH = UNDEF                    # sick from who? 
    nextday_meetcnt::Int16 = 0                  # number of contacts per day
    grade::Int16 = 0                            # school grade 9 - 12
    association::Int16 = 0                      # student association, 1 => general body, 2 => athletics
    mask::Bool = true                           # wearing a mask? 
    tis::Int16   = 0                            # time in state 
    exp::Int16   = 0                            # max statetime
    dur::NTuple{4, Int8} = (0, 0, 0, 0)         # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Int16   = 999                          # day of infection.
    iso::Bool = false                           # isolated (dosn't enter school)
end
Base.show(io::IO, ::MIME"text/plain", z::Student) = dump(z)

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0       
    frelasymp::Float64 = 0.11                   # relative transmission of asymptomatic
    popsize::Int64 = 1000                       # total population size
    calibration::Bool = false                   # calibration mode
    modeltime::Int64 = 14                       # time period
    asymp_props::Float64 = 0.25                 # proportion of 0 - 19 years old that are asymptomatic
end

## load empty constants  (will be loaded for each process if using pmap, monitor for memory issues)
const humans = Array{Student}(undef, 0) 
const p = ModelParameters()  ## setup default parameters


function runsequential() 
    # run n simulations sequentially and average the results. 
    N = 50
    mp = ModelParameters()
    allresults = map(1:N) do x 
        runsim(x, mp)
    end
    sum(allresults) ./ N
end

function runsim(simnum, ip::ModelParameters)
    # runs a single simulation and creates a dataframe
    Random.seed!(simnum*726)
    results = main(ip)            

    # get incidence/prevalnce data
    all = _collectdf(results)
    #insertcols!(all, 1, :sim => simnum)
    return all
end

function main(ip::ModelParameters)
    # reset the parameters for the simulation scenario
    reset_params(ip)
    initialize() # initialize population

    hmatrix = zeros(Int16, p.popsize, p.modeltime)

    # insert initial infected agents into the model
    if p.calibration
        error("calibration code not implemented")
    end
    # split population in agegroups 
    grps = map(x -> findall(y -> y.grade == x, humans), 1:4)

    # start the time loop
    for st = 1:p.modeltime
        @info "time: $st"
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        time_update()
        dyntrans(grps) 
        daily_inf()
        # end of day
    end
    return hmatrix ## return the model state as well as the age groups. 
end

## Initialization Functions 
reset_params() = reset_params(ModelParameters())
function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters, copy the values from ip to p. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end
end
export reset_params

function initialize() 
    # error checks
    p.popsize == 0 && error("no population size given")
    println("popultion size $(p.popsize)")
    resize!(humans, p.popsize) # resize the human array to change population size
    for i = 1:p.popsize 
        humans[i] = Student()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        x.grade = rand(1:4)
        x.exp = 999  ## susceptible (by default) people don't expire.
        x.dur = sample_epi_durations() # sample epi periods   
        # initialize the next day counts (this is important in initialization since dyntrans runs first)
        get_nextday_counts(x)
    end
end

## Data Collection/ Model State functions
function _get_model_state(t, hmatrix)
    # collects the model state at time t
    for i=1:length(humans)
        hmatrix[i, t] = Int(humans[i].health)
    end    
end

function _collectdf(hmatrix)
    ## takes the output of the humans x time matrix and processes it into a dataframe
    #_names_inci = Symbol.(["lat_inc", "mild_inc", "miso_inc", "inf_inc", "iiso_inc", "hos_inc", "icu_inc", "rec_inc", "ded_inc"])    
    #_names_prev = Symbol.(["sus", "lat", "mild", "miso", "inf", "iiso", "hos", "icu", "rec", "ded"])
    mdf_inc, mdf_prev = _get_incidence_and_prev(hmatrix)
    #mdf = hcat(mdf_inc, mdf_prev)    
    #_names_inc = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_INC"))
    #_names_prev = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_PREV"))
    #_names = vcat(_names_inc..., _names_prev...)
    #datf = DataFrame(mdf, _names)
    #insertcols!(datf, 1, :time => 1:p.modeltime) ## add a time column to the resulting dataframe
    return mdf_inc
end

function _get_incidence_and_prev(hmatrix)
    cols = instances(HEALTH)[1:end - 1] ## don't care about the UNDEF health status
    inc = zeros(Int64, p.modeltime, length(cols))
    pre = zeros(Int64, p.modeltime, length(cols))
    for i = 1:length(cols)
        inc[:, i] = _get_column_incidence(hmatrix, cols[i])
        pre[:, i] = _get_column_prevalence(hmatrix, cols[i])
    end
    return inc, pre
end

function _get_column_incidence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for r in eachrow(hmatrix)
        idx = findfirst(x -> x == inth, r)
        if idx !== nothing 
            timevec[idx] += 1
        end
    end
    return timevec
end

function _get_column_prevalence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for (i, c) in enumerate(eachcol(hmatrix))
        idx = findall(x -> x == inth, c)
        if idx !== nothing
            ps = length(c[idx])    
            timevec[i] = ps    
        end
    end
    return timevec
end

function _count_infectors()     
    pre_ctr = asymp_ctr = mild_ctr = inf_ctr = 0
    for x in humans 
        if x.health != SUS ## meaning they got sick at some point
            if x.sickfrom == PRE
                pre_ctr += 1
            elseif x.sickfrom == ASYMP
                asymp_ctr += 1
            elseif x.sickfrom == MILD
                mild_ctr += 1 
            elseif x.sickfrom == INF 
                inf_ctr += 1 
            else 
                error("sickfrom not set right: $(x.sickfrom)")
            end
        end
    end
    return (pre_ctr, asymp_ctr, mild_ctr, inf_ctr)
end

## Disease Functions
function daily_inf()
    # everywhere there is a probability that one of the susceptible individuals will come in as asymptomatic
    dailyprob = 0.005 # move to ModelParameters
    allsus = findall(x -> x.health == SUS, humans)
    for x in allsus
        if rand() < dailyprob
            @info ("debug: human $x will be infected next day")
            move_to_latent(humans[x])
            humans[x].tis = humans[x].exp # force the switch to asymp/pre
        end
    end
end

function sample_epi_durations()
    # when a person is sick, samples the 
    lat_dist = truncated(LogNormal(log(5.2), 0.1), 4, 7) # truncated between 4 and 7
    pre_dist = truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3
    asy_dist = Gamma(5, 1)
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    latents = Int.(round.(rand(lat_dist)))
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    asymps = Int.(ceil.(rand(asy_dist)))
    infs = Int.(ceil.(rand(inf_dist)))
    return (latents, asymps, pres, infs)
end

# Main Time Step of the Model
function time_update()
    # counters to calculate incidence at each time step
    lat=0; pre=0; asymp=0; inf=0; rec=0; ded=0;
    for x in humans 
        x.tis += 1 
        x.doi += 1 # increase day of infection. variable is garbage until person is latent
        if x.tis >= x.exp             
            @match Symbol(x.swap) begin
                :LAT  => begin move_to_latent(x); lat += 1; end
                :PRE  => begin move_to_pre(x); pre += 1; end
                :ASYMP => begin move_to_asymp(x); asymp += 1; end
                :INF  => begin move_to_inf(x); inf +=1; end    
                :REC  => begin move_to_recovered(x); rec += 1; end
                :DED  => begin move_to_dead(x); ded += 1; end
                _    => error("swap expired, but no swap set.")
            end
        end
        # get the meet counts for the next day 
        get_nextday_counts(x)
    end
    return (lat, asymp, inf, rec, ded)
end

function move_to_latent(x::Student)
    ## transfers human h to the incubation period and samples the duration
    x.health = LAT
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 
    x.exp = x.dur[1] # get the latent period
    x.swap = rand() <  p.asymp_props ? ASYMP : PRE 
    ## in calibration mode, latent people never become infectious.
    if p.calibration 
        x.swap = LAT 
        x.exp = 999
    end
end

function move_to_asymp(x::Student)
    ## transfers human h to the asymptomatic stage 
    x.health = ASYMP     
    x.tis = 0 
    x.exp = x.dur[2] # get the presymptomatic period
    x.swap = REC 
end

function move_to_pre(x::Student)
    x.health = PRE
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period
    x.swap = INF
end

function move_to_inf(x::Student)
    ## transfers human h to the severe infection stage for γ days 
    ## simplified function for calibration/general purposes
    x.health = INF
    x.tis = 0 
    x.exp = x.dur[4]
    x.iso = true
    x.swap = REC 
end

function move_to_dead(h::Student)
    # no level of alchemy will bring someone back to life. 
    h.health = DED
    h.swap = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    h.iso = true # a dead person is isolated
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end

function move_to_recovered(h::Student)
    h.health = REC
    h.swap = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely         
end

## Transmission Dynamics 
@inline function _get_betavalue(xhealth) 
    # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
    bf = p.β
    if xhealth == ASYMP
        bf = bf * p.frelasymp
    elseif xhealth == INF 
        bf = bf * 0.89
    end
    return bf
end

@inline function get_nextday_counts(x::Student)
    # gets the number of agents x will meet (or recieve contacts from) 
    # we can sample this at the start of the simulation to avoid everyday but its type stable and fast
    cnt = 0
    ag = x.grade
    x.health != DED && (cnt = rand(nbs[ag]))  # expensive operation, try to optimize
    x.nextday_meetcnt = cnt
    return cnt
end

@inline function beta_mask_adjust(beta, inf_mask, sus_mask) 
    # https://msphere.asm.org/content/5/5/e00637-20
    nb = beta
    if inf_mask && !sus_mask
        nb = beta * rand(27:42) / 100   # 58 to 73% reduction if infected person wearing mask
    elseif !inf_mask && sus_mask 
        nb = beta * rand(63:83) / 100   # 20 - 40% reduction if susceptible person wearing mask
    elseif inf_mask && sus_mask 
        nb = beta * rand(24:29) / 100   # if both of them are wearing mask, figure 2E of the refernece
    end
    return nb
end

@inline function beta_vax_adjust() 
    error("not implemented")
end

function dyntrans(grps)
    totalmet = 0 # count the total number of contacts (total for day, for all INF contacts)
    totalinf = 0 # count number of new infected 
    
    ## find all the people infectious in school. 
    infs = findall(x -> x.health in (PRE, ASYMP), humans)
   
    # go through every infectious person
    for xid in infs 
        x = humans[xid]
        xhealth = x.health
        beta = _get_betavalue(xhealth) 
        cnts = x.nextday_meetcnt        ## how many contacts to give
        cnts == 0 && continue        
        gpw = Int.(round.(cm[x.grade]*cnts)) # split the counts over students in other grades
        for (i, g) in enumerate(gpw)
            # sample the people from each group
            meet = rand(grps[i], g)
            # go through each person
            for j in meet 
                y = humans[j]
                ycnt = y.nextday_meetcnt             
                ycnt == 0 && continue               
                y.nextday_meetcnt = y.nextday_meetcnt - 1 # remove a contact
                totalmet += 1

                # tranmission dynamics
                if  y.health == SUS && y.swap == UNDEF
                    beta = beta_mask_adjust(beta, x.mask, y.mask)
                    if rand() < beta
                        totalinf += 1
                        y.swap = LAT
                        y.exp = y.tis   ## force the move to latent in the next time step.
                        y.sickfrom = xhealth ## stores the infector's status to the infectee's sickfrom                            
                    end  
                end                                    
            end
        end                
    end
    @info("number of infecteds: $(length(infs)), total met $(totalmet), new infections $totalinf")
    return totalmet, totalinf
end
export dyntrans

## Contact Matrix Functions 
function contact_matrix() 
    # regular contacts, just with 5 age groups. 
    #  0-4, 5-19, 20-49, 50-64, 65+
    CM = Array{Array{Float64, 1}, 1}(undef, 5)
    CM[1] = [0.2287, 0.1839, 0.4219, 0.1116+0.0539]
    CM[2] = [0.0276, 0.5964, 0.2878, 0.0591+0.0291]
    CM[3] = [0.0376, 0.1454, 0.6253, 0.1423+0.0494]
    CM[4] = [0.0242, 0.1094, 0.4867, 0.2723+0.1074]
    CM[5] = [0.0207, 0.1083, 0.4071, 0.2193+0.2446]
    return CM
end

function negative_binomials() 
    ## the means/sd here are calculated using _calc_avgag
    means = [10.21, 16.793, 13.7950, 11.2669]
    sd = [7.65, 11.7201, 10.5045, 9.5935]
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms   
end
const nbs = negative_binomials()
const cm = contact_matrix()


## internal functions to do intermediate calculations
function _calc_avgag(lb, hb) 
    ## internal function to calculate the mean/sd of the negative binomials
    ## returns a vector of sampled number of contacts between age group lb to age group hb
    dists = _negative_binomials_15ag()[lb:hb]
    totalcon = Vector{Int64}(undef, 0)
    for d in dists 
        append!(totalcon, rand(d, 10000))
    end    
    return totalcon
end
