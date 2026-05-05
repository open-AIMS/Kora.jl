module UnicodePlotsExt

using OrderedCollections
using UnicodePlots, Term
using Kora

"""
Helper function to generate validation plots for survival models.
"""
function Kora.viz.survival_performance_plots(
    groupings::OrderedDict,
    model_fits::Kora.PolySurvivalModel;
    target_groups::Vector{String}=Kora.TARGET_GROUPS
)
    bin_id_col = Kora.BIN_ID

    for group_id in 1:length(target_groups)
        group_df = groupings[target_groups[group_id]]
        max_group_size = maximum(group_df.diam)

        bin_ids = unique(group_df[!, bin_id_col])

        # Training data plot
        mean_train = [maximum(group_df[group_df[!, bin_id_col].==i, Kora.TRAIN_CLASS_MEAN_ID]) for i in bin_ids]
        std_train = [maximum(group_df[group_df[!, bin_id_col].==i, Kora.TRAIN_CLASS_STD_ID]) for i in bin_ids]

        score_text = build_metric_display(model_fits.performance.train, group_id)
        train_title = "$(target_groups[group_id])\nTraining Data\n$(score_text)"

        plt_train = scatterplot(mean_train; ylim=(0, 1.1), title=train_title)
        lineplot!(plt_train, min.(mean_train .+ std_train, 1.0))
        lineplot!(plt_train, mean_train .- std_train)
        lineplot!(plt_train, [model_fits[group_id](Float64(i)) for i in 0:max_group_size])
        xlabel!(plt_train, "Diameter Bin")
        ylabel!(plt_train, "Survival")

        # Test data plot
        mean_test = [maximum(group_df[group_df[!, bin_id_col].==i, Kora.TEST_CLASS_MEAN_ID]) for i in bin_ids]
        std_test = [maximum(group_df[group_df[!, bin_id_col].==i, Kora.TEST_CLASS_STD_ID]) for i in bin_ids]

        score_text = build_metric_display(model_fits.performance.test, group_id)
        test_title = "$(target_groups[group_id])\nTest Data\n$(score_text)"

        plt_test = scatterplot(mean_test; ylim=(0, 1.1), title=test_title, name="Observed")
        lineplot!(plt_test, min.(mean_test .+ std_test, 1.0), name="+1 stdev")
        lineplot!(plt_test, mean_test .- std_test, name="-1 stdev")
        lineplot!(plt_test, [model_fits[group_id](Float64(i)) for i in 0:max_group_size], name="Model")
        xlabel!(plt_test, "Diameter Bin")
        display(UnicodePlots.panel(plt_train) * UnicodePlots.panel(plt_test))
    end
end

function build_metric_display(results::NamedTuple, idx::Int64)
    res_text = []
    for m in Kora.ALL_METRICS
        _s = getfield(results, Symbol(m))[idx]
        push!(res_text, "$(string(m)): $(round(_s; digits=3))")
    end

    return join(res_text, " | ")
end

"""
Helper function to generate validation plots for growth models.
"""
function Kora.viz.growth_performance_plots(
    groupings::OrderedDict,
    model_fits::Kora.PolyGrowthModel;
    target_groups::Vector{String}=Kora.TARGET_GROUPS
)
    for group_id in 1:length(target_groups)
        group_df = groupings[target_groups[group_id]]

        # Training data plot
        sub_df = group_df[group_df[!, Kora.TRAIN_CLASS].>0, :]
        x_idx = sortperm(sub_df.diam)
        xi = sub_df.diam[x_idx]
        yi = sub_df.diamnext[x_idx]

        score_text = build_metric_display(model_fits.performance.train, group_id)
        train_title = "$(target_groups[group_id])\nTraining Data\n$(score_text)"

        model = model_fits[group_id]
        plt_train = scatterplot(yi; title=train_title)
        lineplot!(plt_train, [model(Float64(i)) for i in xi])
        xlabel!(plt_train, "Diameter [cm]")
        ylabel!(plt_train, "Diameter at t+1 [cm]")

        # As above, but for test data
        sub_df = group_df[group_df[!, Kora.TEST_CLASS].>0, :]
        x_idx = sortperm(sub_df.diam)
        xi = sub_df.diam[x_idx]
        yi = sub_df.diamnext[x_idx]

        score_text = build_metric_display(model_fits.performance.test, group_id)
        test_title = "$(target_groups[group_id])\nTest Data\n$(score_text)"
        plt_test = scatterplot(yi; title=test_title, name="Observed")
        lineplot!(plt_test, [model(Float64(i)) for i in xi], name="Model")
        xlabel!(plt_test, "Diameter [cm]")
        display(UnicodePlots.panel(plt_train) * UnicodePlots.panel(plt_test))
    end
end

end