function p(x, y, X, Y, Ker::Kernels)
    A = [collect(x);collect(y)]
    isempty(Y) && return kernel(A, [i[1:length(x)] for i in X], Ker)
    B = [[i[1:length(x)];j[1:length(y)]] for (i,j) in zip(X,Y)]
    #I give to the kernel A which is the vector with the X and Y value of timeseries i want to know the probability of (1 x k + l)
    #and B wich is the vector of vectors of all realization of time series
    return kernel(A, B, Ker)
end

function kernel(A, B, x::Frequency)
    count = 0
    for i in B
        i == A && (count += 1)
    end
    #println(count)
    return count/length(B)
end

function Time_series(Dynamic, n, Time)
    #Dynamic is the struct that define the time series
    #n is the number of time steps
    if length(Dynamic.x) == Time
        x = Dynamic.x
    else
        x = [i for j in 1:(Time÷length(Dynamic.x))+1 for i in Dynamic.x]
    end
    if length(Dynamic.y) == Time
        y = Dynamic.y
    else
        y = [i for j in 1:Time÷length(Dynamic.y)+1 for i in Dynamic.y]
    end


    X = Vector{Vector{Float64}}(undef,n)
    Y = Vector{Vector{Float64}}(undef,n)
    for j in 1:n
        X[j] = zeros(Float64, Time)
        Y[j] = zeros(Float64, Time)
        X[j][1] = Dynamic.x0(Dynamic.x)
        Y[j][1] = Dynamic.y0(Dynamic.y)
        for i in 2:Time
            X[j][i] = Dynamic.f(X[j][1:i-1], Y[j][1:i-1], Dynamic.a...)
            Y[j][i] = Dynamic.g(X[j][1:i-1], Y[j][1:i-1], Dynamic.b...)
        end
    end
    return Time_Series(n, x, y, X, Y, Dynamic)
end

#inizialization function

function init_rand(c)
    return rand(c[1])
end

function full_prob(t,k,l,Time_series, Ker::Kernels)
    X = collect(Iterators.product(Time_series.x[t-k:t]...))
    Y = collect(Iterators.product(Time_series.y[t-l:t-1]...))
    I = [i[t-k:t] for i in Time_series.X]
    J = [j[t-l:t-1] for j in Time_series.Y]
    B = [[i;j] for (i,j) in zip(I,J)]
    H = [collect((x...,y...)) for x in X, y in Y]
    return [kernel(i, B, Ker) for i in H]
end

function full_prob_k(full, K, k, L, l,xdim,ydim)
    kdiff = K-k
    ldiff = L-l
    Z = collect(Iterators.product([1:xdim for i in 1:kdiff]...))
    X = collect(Iterators.product([1:ydim for i in 1:ldiff]...))
    kcycle = [collect(z) for z in Z]
    lcycle = [collect(x) for x in X]
    kfull = [1:xdim for i in 1:k+1]
    Lfull = [1:ydim for i in 1:L]
    lfull = [1:ydim for i in 1:l]
    if kdiff == 0
        A = full
    else
        A = sum([full[i...,kfull...,Lfull...] for i in kcycle])
    end
    if ldiff == 0
        return A[kfull...,lfull...]
    else
        B = sum([A[kfull...,i...,lfull...] for i in lcycle])
    end
    return B[kfull...,lfull...]
end

function marginal_x_t_y_t(full, K, L,xdim,ydim)
    Kfull = [1:xdim for i in 1:K]
    Lfull = [1:ydim for i in 1:L]
    return sum([full[Kfull...,i,Lfull...] for i in 1:xdim])
end

function TE_dist(full, K, L,xdim,ydim)
    full_px = full_prob_k(full, K, K, L, 0,xdim,ydim)
    pxy = marginal_x_t_y_t(full, K, L,xdim,ydim)
    px = full_prob_k(full, K-1, K-1, L+1, 0,xdim,ydim)
    M = zeros(Float64, length(full), 4)
    M[:,1] = full[1:end]
    M[:,2] = repeat(full_px[1:end], outer = Int64(length(full)/length(full_px)))
    M[:,4] = repeat(px[1:end], outer = Int64(length(full)/length(px)))
    M[:,3] = reduce(vcat, [repeat(pxy[(k-1)*length(px)+1:(k)*length(px)], outer = 2) for k in 1:Int64(length(pxy)/length(px))])
    return M[M[:,1].!=0,:]
end

function TE_value(Dist)
    sum = 0.0
    for i in 1:size(Dist)[1]
        sum += Dist[i,1]*log2(Dist[i,1]*Dist[i,4]/(Dist[i,2]*Dist[i,3]))
    end
    return sum
end


function TE(t, k, l ,Time_series::Time_Series, Ker::Kernels)
    K = k[end]
    L = l[end]
    xlen = length(Time_series.x[1])
    ylen = length(Time_series.y[1])
    Dist = full_prob(t,K,L,Time_series, Ker)
    V = []
    for i in k
        for j in l
            A = TE_dist(full_prob_k(Dist, K, i, L, j,xlen,ylen),i,j,xlen,ylen)
            push!(V, TE(i, j, Time_series.n, TE_value(A), t, Time_series.Dynamic, Ker, A))
        end
    end
    return V
end

#Various Dynamical functions
copy() = println("coefficent of this function should be l. return the value of y_t-l")
copy(x,y,l) = length(y) > L && (return y[end-l])

function both_sum()
    println("coefficent of this function should be, k, l, xval. If the sum of x_(t-k) and y_(t-l) is greater than 1 return xval[1] else xval[end]")
end

function both_sum(x, y, k , l , xval)
    length(y) < l+1 && return 0
    length(x) < k+1 && return 0
    x[end-k]+y[end-l] >= 1 ? (return xval[1]) : (return xval[end])
end

function rand_det()
    println("coefficent of this function should be, p, l, Xval, Coeff. Coeff length should be equal to l")
end

function rand_det(x,y,p,l,Xval, Coeff)
    length(y) < l+1 && return 0
    if rand()  < p
        return rand(Xval)
    else
        Det_evolution(x, y, l, Coeff...)
    end
end

function rand_on_value()
    println("coefficent of this function should be, p, xval; return random value with probability p, else copy the last y")
end
    
function rand_on_value(x,y,p,xval)
    length(y) < 2 && return 0
    rand() < p ? (return rand(xval)) : (return y[end-1])
end

function Double_dependence()
    println("coefficent of this function should be, k_1, l_1, k_2, l_2, xval. If both both sum with different l and k give \n the same value return xval[1] else xval[end]")
end

function Double_dependence(x, y, k_1 ,l_1 ,k_2,l_2, xval)
    both_sum(x, y, k_1,l_1, xval) == xval[1] && both_sum(x, y, k_2,l_2, xval) == xval[1] ? (return xval[1]) : (return xval[end])
end

function rand_copy()
    println("coefficent of this function should be, p, xval: probability of copying value of last y, 1-p it flips it")
end

function rand_copy(x,y,p, xval)
    length(y) < 2 && return rand(xval)
    rand() < p ? (return y[end-1]) : (return -y[end-1])
end

function Det_evolution()
    println("coefficent of this function should be, l, Coeff. Coeff length should be equal to l. if the sum of the last l values of y 
    \n multiplied by the coefficients is >= 1 then it returns 1, else -1")
end

function Det_evolution(x, y, l,Coeff...)
    length(y) < l && return 1
    C = [i for i in Coeff]
    sum(y[end-l+1:end].*C)>=1 ? (return 1) : (return -1)
end

ra() = println("coefficent of this function should be c, the values it can takes randomly")
ra(x,y, c...) = rand(c)

flip(x,y, c...) = -y[end];

xor() = println("coefficent of this function should be xval. If the last value of x and y are the same return xval[1] else xval[end]")

function xor(x,y,xval)
    length(y) < 2 && return 0
    (x[end] == y[end]) ? (return xval[1]) : (return xval[end])
end

double_series() = println("coefficent of this function should be xval. If the last value of y and the series Z (it's a random series \n 
so just write rand(xval)) are the same return xval[1] else xval[end]")

function double_series(x,y,xval)
    length(y) < 2 && return rand(xval)
    (rand(xval) == y[end]) ? (return xval[1]) : (return xval[end])
end

p_rand_yev() = println("coefficent of this function should be p, yval. If the probability p is greater than a random number return \n
a random value from yval else return the last value of y multiplied by -1")

function p_rand_yev(x,y,p,yval)
    length(y) < 2 && return rand(yval)
    rand() < p ? (return rand(yval)) : (return -y[end-1])
end

