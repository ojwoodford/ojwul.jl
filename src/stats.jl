export samplenpdistribution

function samplenpdistribution(array, num=1)
    cdf = cumsum(vec(array))
    ind = searchsortedlast.(Ref(cdf), rand(num) * cdf[end])
    if ndims(array) > 1
        i2s = CartesianIndices(size(array))
        return @inbounds i2s[ind]
    else
        return ind
    end
end

