using Images

# Function definitions
"Load a set of similarly-named images"
function get_imageset(dir, prefix)
    return [load(joinpath(dir, file)) for file in readdir(dir) if !isnothing(match(Regex("$(prefix).*00.*\\.tiff\$"), file))]
end
"Compute the distance between the bottom of the first image and the top of the second image, 
given an overlap"
function overlap_dist(im1, im2, guess)
    footer = size(im2, 1) ÷ 17
    # display(hcat(im1[imsize-guess:imsize, :], im2[1:guess+1, :]))
    mae((@view im1[end-footer-guess:end-footer, :]), (@view im2[1:guess+1, :]))
end
"""
Identify the best overlap between two images, then stich together at the halfway point.
Recursively applied for more than two images.
"""
function vert_stitch(args...; marking=false)
    if length(args) == 1
        return args[1]
    elseif length(args) == 2
        # return vert_stitch2(args[1], args[2]; marking)
        im1, im2 = args
        footer = size(im2, 2) ÷ 17
        perc05 = size(im2, 1) ÷ 20
        possible_overlaps = perc05:(size(im2, 1)-footer-3perc05)
        best_overlap = possible_overlaps[argmin([overlap_dist(im1, im2, guess) for guess in possible_overlaps])]
        half1_overlap = ÷(best_overlap, 2)
        half2_overlap = best_overlap - half1_overlap
        @info "stitch" size(im1) size(im2) best_overlap overlap_dist(im1, im2, best_overlap)
        stitched = vcat(
            im1[1:end-half1_overlap-footer, :],
            fill(Gray{N0f8}(1), (marking, size(im1, 2))), # Put a white line at stitch, if desired
            im2[half2_overlap:end, :],
        )
        return stitched
    else
        stitched = args[1]
        for next in args[2:end]
            stitched = vert_stitch(stitched, next; marking)
        end
        return stitched
    end
end

# Identify data directories and prefixes
dirs = ["Fulbright_SEM_batch3", "Fulbright_SEM_batch4", "Fulbright_SEM_batch5"]
prefixes = Dict(
    dirs[1] => ["1_center_top", "1_right_top", "2_center_top", "2_left_top", "3_center_top", "3_left_top", "4_center_top", "4_right_top", "5_center_top", "5_left2_top", "6_center_top", "6_left_top"],
    dirs[2] => ["1_top_center", "1_top_right", "2_top_center", "2_top_right", "6_top_right", "6_top_left"],
    dirs[3] => ["1_mid_top", "2_left", "2_mid_top", "4_mid", "5_mid", "5_right", "6_mid", "6_right"],
)
im1ex = load(datadir("exp_raw", dirs[1], readdir(datadir("exp_raw", dirs[1]))[3]))
im2ex = load(datadir("exp_raw", dirs[1], readdir(datadir("exp_raw", dirs[1]))[4]))
@time mae(im1ex, im2ex)
@time overlap_dist(im1ex, im2ex, 100)

for batch in dirs
    imsets = prefixes[batch]
    for imset in imsets
        fname = "$(batch)_$(imset)_stitched.tiff"
        if isfile(datadir("exp_raw", "merged_SEM"), fname)
            @info "Skipping $(batch) with prefix $(imset), already stitched."
            continue
        end
        set = get_imageset(datadir("exp_raw", batch), imset)
        @info "Stitching images in $(batch) with prefix $(imset)" size(set)
        @time stitched = vert_stitch(set...; marking=false)
        save(datadir("exp_raw", batch, "$(imset)_stitched_image.tiff"), stitched)
        save(datadir("exp_raw", "merged_SEM", fname), stitched)
    end
end
# Special cases with manual fixups:
# batch 3, 5_center_top
# batch 5, 2_mid_top


# Listing of experimental conditions
conditions_2 = Table((
    batch=[2, 2, 2, 2 ],
    semno=[1, 2, 3, 4],
    batchno=[12, 7, 18, 2],
    rundate=["2025-02-25", "2025-02-25", "2025-02-27", "2025-02-27"],
    formulation = ["sucrose 10%", "sucrose 10%", "sucrose 10%", "sucrose 10%"],
    Tsh = [-20, -20, -30, -30]u"°C",
))
conditions_3 = Table((
    batch=[3, 3, 3, 3, 3, 3 ],
    semno=[1, 2, 3, 4, 5, 6],
    batchno=[8, 10, 14, 18, 2, 3],
    rundate=["2025-01-20", "2025-01-20", "2025-01-28", "2025-01-28", "2025-02-03", "2025-02-03"],
    formulation = ["sucrose 5%", "sucrose 5%", "sucrose 5%", "sucrose 5%", "sucrose 5%", "sucrose 5%", ],
    Tsh = [-20, -20, -10, -10, 0, 0]u"°C",
))
conditions_4 = Table((
    batch = [4, 4, 4],
    semno = [ 1, 2, 6,],
    batchno=[7, 3, 19],
    rundate = ["2025-05-20", "2025-05-20", "2025-05-23"],
    formulation=["sucrose 5% CIN", "sucrose 5% CIN", "mannitol 5%", ],
    Tsh = [-20, -20, 10]u"°C",
))
conditions_5 = Table((
    batch = [5, 5, 5, 5, 5],
    semno = [1, 2, 4, 5, 6],
    batchno=[14, 18, 4, 10, 13],
    rundate = ["2025-06-09", "2025-06-09", "2025-06-16", "2025-06-18", "2025-06-18"],
    formulation = ["dextran 5% 2mL", "dextran 5% 2mL", "dextran 5% 3mL", "dextran 5% 3mL", "dextran 5% 3mL"],
    Tsh = [5, 5, -10, 0, 0]u"°C",
))

conditions = vcat(conditions_2, conditions_3, conditions_4, conditions_5)

using CSV
using CSVFiles
safesave(datadir("exp_raw", "merged_SEM", "lyo_conditions.csv"),  conditions)