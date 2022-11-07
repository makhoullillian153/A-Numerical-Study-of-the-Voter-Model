#import Pkg
#Pkg.add("StatsBase")
using Distributions
using StatsBase

function neighbor(nbhd,i,j,boundaries = true)
    
    a = 0
    b = 0

    # possible movements
    movements = [0 1; 0 -1 ;-1 0; 1 0]

    if boundaries
        while a <= 0 || b <= 0
            index = sample(1:1:4,1)[1]
            a = movements[index,1]
            b = movements[index,2]

            if a > size(nbhd)[1] || a > size(nbhd)[2]
                a = 0
                b = 0
            end
        end        
    else
        index = sample(1:1:4,1)
        a = movements[1,index]
        b = movements[2,index]
        
        if a == 0
            a = size(nbhd)[1]
        end
        if a == size(nbhd)[1] + 1
            a = 1
        end
        if b == 0
            b = size(nbhd)[2]
        end
        if b == size(nbhd)[2] + 1
            b = 1
        end
    end

    return [a,b]
end

function V4(row, col, s, pOne = 0.5)
    
    N = row*col
    consensusT = Vector{Int128}
    for k in 1:s
        states = sample([0,1], weights([1-pOne,pOne]) , N)
        nbhd = reshape(states, row,col)

        while true
            pChange = [N-sum(nbhd)/N, sum(nbhd)/N]

            if sum(nbhd) == 0 || sum(nbhd) == N
                break
            end

            # randomly select an individual
            i = sample(1:1:row,1)[1]
            j = sample(1:1:col,1)[1]

            # select a neighbor
            lst = neighbor(nbhd,i,j)
            a = lst[1]
            b = lst[2]

            if pChange[nbhd[i,j] + 1] > rand(1)
                nbhd[i,j] = nbhd[a,b]
            end

            consensusT[k] += 1
        end
    end
    return consensusT
end

println(mean(V4(2,2,1000)))