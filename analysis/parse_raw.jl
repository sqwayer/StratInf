using FileIO, JSON, DataFrames, CSV, ProgressMeter

function parse_raw(foldername)
    filesList = filter(x -> occursin(".txt", x), readdir(string("WMM/", foldername)))

    wb = Progress(length(filesList), 1)
    sessPerSub = Dict()
    for fname in filesList

        open(string("WMM/", foldername,"/", fname), "r") do f
            D = JSON.parse(f)
        end
        D = filter(di -> haskey(di, "part") && di["part"] == "stim" && !di["training"], D)
        df = DataFrame(merge(vcat, D...))
        for n in names(df)
            if any(isnothing.(df[!,n]))
                df[!,n] = convert(Vector{Union{eltype(df[!,n]), Missing}}, df[!,n])
                replace!(df[!,n], nothing => missing)
            end
        end
        # Correct session numbers
        if haskey(sessPerSub, df.subject[1])
            sessNum = sessPerSub[df.subject[1]] + 1
            sessPerSub[df.subject[1]] += 1
        else
            sessNum = 0
            sessPerSub[df.subject[1]] = 0
        end
        csvname = string("WMM/", foldername,"/", df.subject[1], "_", sessNum, ".csv")
        CSV.write(csvname, df)
        next!(wb)
    end

    return
end
