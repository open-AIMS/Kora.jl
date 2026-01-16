using CoralFlow.Interpolations

"""
Helper function to generate validation plots for survival models using Makie.
"""
function CoralFlow.viz.survival_performance_plots(
    groupings::OrderedDict,
    model_fits::CoralFlow.PolySurvivalModel;
    target_groups::Vector{String}=CoralFlow.TARGET_GROUPS,
    figsize=(1200, 400),
    save_path=nothing
)
    bin_id_col = CoralFlow.BIN_ID

    for group_id in eachindex(target_groups)
        group_df = groupings[target_groups[group_id]]
        min_group_size = minimum(group_df.diam)
        max_group_size = maximum(group_df.diam)

        bin_ids = sort(unique(group_df[!, bin_id_col]))

        # Training data
        mean_train = [
            maximum(group_df[group_df[!, bin_id_col] .== i, CoralFlow.TRAIN_CLASS_MEAN_ID])
            for i in bin_ids
        ]
        std_train = [
            maximum(group_df[group_df[!, bin_id_col] .== i, CoralFlow.TRAIN_CLASS_STD_ID])
            for i in bin_ids
        ]

        # Test data
        mean_test = [
            maximum(group_df[group_df[!, bin_id_col] .== i, CoralFlow.TEST_CLASS_MEAN_ID])
            for i in bin_ids
        ]
        std_test = [
            maximum(group_df[group_df[!, bin_id_col] .== i, CoralFlow.TEST_CLASS_STD_ID])
            for i in bin_ids
        ]

        # Map bin IDs to actual diameter values and ranges for adaptive bins
        bin_centers = Float64[]
        bin_mins = Float64[]
        bin_maxs = Float64[]

        for bin_id in bin_ids
            bin_data = group_df[group_df[!, bin_id_col] .== bin_id, :diam]
            push!(bin_centers, (minimum(bin_data) + maximum(bin_data)) / 2)
            push!(bin_mins, minimum(bin_data))
            push!(bin_maxs, maximum(bin_data))
        end

        # Model predictions across full diameter range
        model_x = range(min_group_size, max_group_size; length=100)
        model_y = [model_fits[group_id](Float64(x)) for x in model_x]

        # Create interpolation functions
        mean_train_interp = linear_interpolation(
            bin_centers, mean_train; extrapolation_bc=Flat()
        )
        std_train_interp = linear_interpolation(
            bin_centers, std_train; extrapolation_bc=Flat()
        )
        mean_test_interp = linear_interpolation(
            bin_centers, mean_test; extrapolation_bc=Flat()
        )
        std_test_interp = linear_interpolation(
            bin_centers, std_test; extrapolation_bc=Flat()
        )

        # Calculate interpolated values
        mean_train_full = mean_train_interp.(model_x)
        std_train_full = std_train_interp.(model_x)
        mean_test_full = mean_test_interp.(model_x)
        std_test_full = std_test_interp.(model_x)

        # Create figure with side-by-side subplots and space for legend
        fig = Figure(; size=figsize)

        # Plot data for both training and test
        plot_data = [
            (
                mean_train,
                std_train_full,
                mean_train_full,
                "Training Data",
                model_fits.performance.train
            ),
            (
                mean_test,
                std_test_full,
                mean_test_full,
                "Test Data",
                model_fits.performance.test
            )
        ]

        plot_handles = []

        group_name = replace(titlecase(target_groups[group_id]), "_" => " ")
        for (col_idx, (observed_means, std_full, mean_full, data_type, performance)) in
            enumerate(plot_data)
            ax = Axis(fig[1, col_idx];
                title="$(group_name) - $data_type\n$(build_metric_display(performance, group_id))",
                xlabel="Diameter [cm]",
                ylabel="Survival",
                limits=(nothing, nothing, 0, 1.1)
            )

            # Plot observed data with error bands and bin range indicators
            p_obs = scatter!(ax, bin_centers, observed_means; color=:blue, markersize=8)

            # Add horizontal error bars to show diameter range of each bin
            p_range = errorbars!(ax, bin_centers, observed_means,
                bin_centers .- bin_mins, bin_maxs .- bin_centers;
                direction=:x, color=(:darkblue, 0.3), linewidth=2)

            # Survival probability error bands
            p_band = band!(ax, model_x,
                max.(mean_full .- std_full, 0),
                min.(mean_full .+ std_full, 1.0);
                color=(:blue, 0.3))

            p_model = lines!(ax, model_x, model_y; color=(:red, 0.5), linewidth=2)

            # Store handles from first plot for legend
            if col_idx == 1
                plot_handles = [p_obs, p_range, p_band, p_model]
            end
        end

        # Add shared legend to the right of the figure
        Legend(fig[1, 3], plot_handles, ["Bin Mean", "Bin Range", "±1 stdev", "Model"])

        _display_or_save(fig, "survival", target_groups[group_id], save_path)
    end
end

function build_metric_display(results::NamedTuple, idx::Int64)
    res_text = []
    for m in CoralFlow.ALL_METRICS
        _s = getfield(results, Symbol(m))[idx]
        push!(res_text, "$(string(m)): $(round(_s; digits=3))")
    end

    return join(res_text, " | ")
end

"""
Helper function to generate validation plots for growth models using Makie.
"""
function CoralFlow.viz.growth_performance_plots(
    groupings::OrderedDict,
    model_fits::CoralFlow.PolyGrowthModel;
    target_groups::Vector{String}=CoralFlow.TARGET_GROUPS,
    figsize=(1200, 400),
    save_path=nothing,
    alpha=0.6
)
    for group_id in eachindex(target_groups)
        group_df = groupings[target_groups[group_id]]

        # Training data
        train_df = group_df[group_df[!, CoralFlow.TRAIN_CLASS] .> 0, :]
        train_x_idx = sortperm(train_df.diam)
        train_xi = train_df.diam[train_x_idx]
        train_yi = train_df.diamnext[train_x_idx]

        # Test data
        test_df = group_df[group_df[!, CoralFlow.TEST_CLASS] .> 0, :]
        test_x_idx = sortperm(test_df.diam)
        test_xi = test_df.diam[test_x_idx]
        test_yi = test_df.diamnext[test_x_idx]

        # Model predictions
        model = model_fits[group_id]
        all_x = vcat(train_xi, test_xi)
        x_range = range(minimum(all_x), maximum(all_x); length=100)
        model_y = [model(Float64(x)) for x in x_range]

        # Create figure with side-by-side subplots and space for legend
        fig = Figure(; size=figsize)

        # Training plot
        group_name = replace(titlecase(target_groups[group_id]), "_" => " ")
        ax1 = Axis(fig[1, 1];
            title="$(group_name) - Training Data\n$(build_metric_display(model_fits.performance.train, group_id))",
            xlabel="Diameter [cm]",
            ylabel="Diameter at t+1 [cm]"
        )

        p1_obs = scatter!(ax1, train_xi, train_yi;
            color=(:blue, alpha), markersize=6)
        p1_model = lines!(ax1, x_range, model_y;
            color=(:red, 0.5), linewidth=2)

        # Test plot
        ax2 = Axis(fig[1, 2];
            title="$(group_name) - Test Data\n$(build_metric_display(model_fits.performance.test, group_id))",
            xlabel="Diameter [cm]",
            ylabel="Diameter at t+1 [cm]"
        )

        p2_obs = scatter!(ax2, test_xi, test_yi;
            color=(:blue, alpha), markersize=6)
        p2_model = lines!(ax2, x_range, model_y;
            color=(:red, 0.5), linewidth=2)

        # Add reference line (y = x) for comparison
        min_val = min(minimum(all_x), minimum(vcat(train_yi, test_yi)))
        max_val = max(maximum(all_x), maximum(vcat(train_yi, test_yi)))
        p1_ref = lines!(ax1, [min_val, max_val], [min_val, max_val];
            color=:gray, linestyle=:dash, alpha=0.5)
        p2_ref = lines!(ax2, [min_val, max_val], [min_val, max_val];
            color=:gray, linestyle=:dash, alpha=0.5)

        # Add shared legend to the right of the figure
        Legend(fig[1, 3], [p1_obs, p1_model, p1_ref], ["Observed", "Model", "y = x"])

        # Display or save
        _display_or_save(fig, "growth", target_groups[group_id], save_path)
    end
end

# Additional utility functions for enhanced Makie plotting

"""
Create a comprehensive dashboard showing both survival and growth models.
"""
function CoralFlow.viz.model_dashboard(
    groupings::OrderedDict,
    survival_fits::CoralFlow.PolySurvivalModel,
    growth_fits::CoralFlow.PolyGrowthModel;
    target_groups::Vector{String}=CoralFlow.TARGET_GROUPS,
    figsize=(1200, 800),
    save_path=nothing
)
    fig = Figure(; size=figsize)

    for (i, group) in enumerate(target_groups)
        # Create survival subplot
        ax_surv = Axis(fig[i, 1];
            title="$(group) - Survival",
            xlabel="Diameter Bin",
            ylabel="Survival Probability"
        )

        # Create growth subplot
        ax_growth = Axis(fig[i, 2];
            title="$(group) - Growth",
            xlabel="Diameter [cm]",
            ylabel="Diameter at t+1 [cm]"
        )

        # Plot survival data (simplified version)
        group_df = groupings[group]
        bin_ids = sort(unique(group_df[!, CoralFlow.BIN_ID]))
        mean_test = [
            maximum(
                group_df[group_df[!, CoralFlow.BIN_ID] .== j, CoralFlow.TEST_CLASS_MEAN_ID]
            ) for j in bin_ids
        ]

        scatter!(ax_surv, bin_ids, mean_test; color=:blue, markersize=6)

        max_group_size = maximum(group_df.diam)
        model_x = bin_ids
        model_y = [survival_fits[i](Float64(j)) for j in model_x]
        lines!(ax_surv, model_x, model_y; color=:red, linewidth=2)

        # Plot growth data (simplified version)
        test_df = group_df[group_df[!, CoralFlow.TEST_CLASS] .> 0, :]
        if !isempty(test_df)
            scatter!(ax_growth, test_df.diam, test_df.diamnext;
                color=(:blue, 0.5), markersize=4)

            x_range = range(minimum(test_df.diam), maximum(test_df.diam); length=50)
            growth_model_y = [growth_fits[i](Float64(x)) for x in x_range]
            lines!(ax_growth, x_range, growth_model_y; color=:red, linewidth=2)
        end
    end

    return _display_or_save(fig, "all", target_groups[group_id], save_path)
end
