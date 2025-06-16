module UnicodePlotsExt

using OrderedCollections
using UnicodePlots, Term
using CoralFlow

"""
Helper function to generate validation plots for survival models.
"""
function CoralFlow.viz.survival_performance_plots(
    groupings::OrderedDict,
    model_fits::CoralFlow.PolySurvivalModel;
    target_groups::Vector{String}=CoralFlow.TARGET_GROUPS
)
    logclass = CoralFlow.LOGCLASS_ID

    for group_id in 1:length(target_groups)
        group_df = groupings[target_groups[group_id]]
        max_group_size = maximum(group_df.diam)

        n_bins = length(unique(group_df[!, logclass]))

        # Training data plot
        mean_train = [maximum(group_df[group_df[!, logclass].==i, CoralFlow.TRAIN_CLASS_MEAN_ID]) for i in 1:n_bins]
        std_train = [maximum(group_df[group_df[!, logclass].==i, CoralFlow.TRAIN_CLASS_STD_ID]) for i in 1:n_bins]

        rmse = round(model_fits.performance.train.rmse[group_id]; digits=3)
        r2 = round(model_fits.performance.train.r2[group_id]; digits=3)
        train_title = "$(target_groups[group_id])\nTraining Data\nRMSE: $(rmse) | R²: $(r2)"

        plt_train = scatterplot(mean_train; ylim=(0, 1.1), title=train_title)
        lineplot!(plt_train, min.(mean_train .+ std_train, 1.0))
        lineplot!(plt_train, mean_train .- std_train)
        lineplot!(plt_train, [model_fits[group_id](Float64(i)) for i in 0:max_group_size])
        xlabel!(plt_train, "Diameter [cm]")
        ylabel!(plt_train, "Survival")

        # Test data plot
        mean_test = [maximum(group_df[group_df[!, logclass].==i, CoralFlow.TEST_CLASS_MEAN_ID]) for i in 1:n_bins]
        std_test = [maximum(group_df[group_df[!, logclass].==i, CoralFlow.TEST_CLASS_STD_ID]) for i in 1:n_bins]

        rmse = round(model_fits.performance.test.rmse[group_id]; digits=3)
        r2 = round(model_fits.performance.test.r2[group_id]; digits=3)
        test_title = "$(target_groups[group_id])\nTest Data\nRMSE: $(rmse) | R²: $(r2)"
        plt_test = scatterplot(mean_test; ylim=(0, 1.1), title=test_title, name="Training data")
        lineplot!(plt_test, min.(mean_test .+ std_test, 1.0), name="+1 stdev")
        lineplot!(plt_test, mean_test .- std_test, name="-1 stdev")
        lineplot!(plt_test, [model_fits[group_id](Float64(i)) for i in 0:max_group_size], name="Model")
        xlabel!(plt_test, "Diameter [cm]")
        display(UnicodePlots.panel(plt_train) * UnicodePlots.panel(plt_test))
    end
end

"""
Helper function to generate validation plots for growth models.
"""
function CoralFlow.viz.growth_performance_plots(
    groupings::OrderedDict,
    model_fits::CoralFlow.PolyGrowthModel;
    target_groups::Vector{String}=CoralFlow.TARGET_GROUPS
)
    for group_id in 1:length(target_groups)
        group_df = groupings[target_groups[group_id]]

        # Training data plot
        sub_df = group_df[group_df[!, CoralFlow.TRAIN_CLASS].>0, :]
        x_idx = sortperm(sub_df.diam)
        xi = sub_df.diam[x_idx]
        yi = sub_df.diamnext[x_idx]

        rmse = round(model_fits.performance.train.rmse[group_id]; digits=3)
        r2 = round(model_fits.performance.train.r2[group_id]; digits=3)
        train_title = "$(target_groups[group_id])\nTraining Data\nRMSE: $(rmse) | R²: $(r2)"

        # Main.@infiltrate
        model = model_fits[group_id]
        plt_train = scatterplot(yi; title=train_title)
        lineplot!(plt_train, [model(Float64(i)) for i in xi])
        xlabel!(plt_train, "Diameter [cm]")
        ylabel!(plt_train, "Diameter at t+1 [cm]")

        # As above, but for test data
        sub_df = group_df[group_df[!, CoralFlow.TEST_CLASS].>0, :]
        x_idx = sortperm(sub_df.diam)
        xi = sub_df.diam[x_idx]
        yi = sub_df.diamnext[x_idx]

        rmse = round(model_fits.performance.test.rmse[group_id]; digits=3)
        r2 = round(model_fits.performance.test.r2[group_id]; digits=3)
        test_title = "$(target_groups[group_id])\nTest Data\nRMSE: $(rmse) | R²: $(r2)"
        plt_test = scatterplot(yi; title=test_title, name="Training data")
        lineplot!(plt_test, [model(Float64(i)) for i in xi], name="Model")
        xlabel!(plt_test, "Diameter [cm]")
        display(UnicodePlots.panel(plt_train) * UnicodePlots.panel(plt_test))
    end
end

end