function genParsGrid(parsRange, parsNet)
# Generate the grid of network parameters.
# Each element in the grid is a Dictionary storing a combination of
# network parameteres
# Wen-Hao Zhang
# Jan-18, 2018

    namePar_MltVal = collect(keys(parsRange))

    sizeGrid = zeros(Int, length(parsRange));
    iter = 1
    for namePar in namePar_MltVal
        sizeGrid[iter] = length(parsRange[namePar])
        iter += 1
    end

    parsGrid = [copy(parsNet) for i = 1:prod(sizeGrid)]
    parsGrid = reshape(parsGrid, tuple(sizeGrid...))


    # Allocate values to the array of composite type
    for iter = 1: length(parsGrid)
        subs = ind2sub(parsGrid, iter);

        iter_name = 1;
        for namePar in namePar_MltVal
            parsGrid[iter][namePar] = parsRange[namePar][subs[iter_name]]
            
            iter_name += 1
        end
    end

    return parsGrid

end
