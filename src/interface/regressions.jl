function fit_growth_models(groupings::OrderedDict{String,DataFrame}; degree::Int=2)::PolyGrowthModel
    models = PolyGrowthFunction[]
    train_rmse = Float32[]
    train_r2 = Float32[]

    test_rmse = Float32[]
    test_r2 = Float32[]

    for (_, df) in groupings
        isempty(df) && continue

        sub_df = df[df[!, TRAIN_CLASS].>0, :]

        x_idx = sortperm(sub_df.diam)
        xi = sub_df.diam[x_idx]
        yi = sub_df.diamnext[x_idx]

        # Fit model
        m = curve_fit(Polynomial, log.(xi), log.(yi), degree)

        # Create growth function
        model = PolyGrowthFunction(xi, yi, m)
        prediction = model.(xi)

        push!(models, model)
        push!(train_rmse, RMSE(prediction, yi))
        push!(train_r2, R2(prediction, yi))

        # As above, but for test data
        sub_df = df[df[!, TEST_CLASS].>0, :]
        x_idx = sortperm(sub_df.diam)
        xi = sub_df.diam[x_idx]
        yi = sub_df.diamnext[x_idx]

        prediction = model.(xi)
        push!(test_rmse, RMSE(prediction, yi))
        push!(test_r2, R2(prediction, yi))
    end

    return PolyGrowthModel(
        collect(keys(groupings)),
        models,
        (;
            train=(; rmse=train_rmse, r2=train_r2),
            test=(; rmse=test_rmse, r2=test_r2),
        )
    )
end

"""
    fit_survival_models(
        groupings::OrderedDict{String, DataFrame}; degree=2
    )

Fit survival models to grouped coral data using logistic regression.

# Returns
Tuple of (models, mcfadden_r2_scores, log_likelihood_scores, brier_scores)
"""
function fit_survival_models(groupings::OrderedDict{String,DataFrame}; degree::Int64=2)::PolySurvivalModel
    models = PolySurvivalFunction[]
    train_rmse = Float32[]
    train_r2 = Float32[]

    test_rmse = Float32[]
    test_r2 = Float32[]
    for (_, df) in groupings
        isempty(df) && continue

        sub_df = df[df.class_train.>0, :]

        x_idx = sortperm(sub_df.diam)
        xi = sub_df.diam[x_idx]
        yi = sub_df.class_train_mean[x_idx]

        # Fit model
        m = curve_fit(Polynomial, log.(xi), yi, degree)

        # Create survival function
        model = PolySurvivalFunction(xi, yi, m)

        # Calculate metrics
        predictions = model.(xi)
        push!(models, model)

        # Save performance scores
        push!(train_rmse, RMSE(predictions, yi))
        push!(train_r2, R2(predictions, yi))

        # Repeat above for test data
        sub_df = df[df.class_test.>0, :]
        x_idx = sortperm(sub_df.diam)
        xi = sub_df.diam[x_idx]
        yi = sub_df.class_test_mean[x_idx]
        predictions = model.(xi)
        push!(test_rmse, RMSE(predictions, yi))
        push!(test_r2, R2(predictions, yi))
    end

    return PolySurvivalModel(
        collect(keys(groupings)),
        models,
        (;
            train=(; rmse=train_rmse, r2=train_r2),
            test=(; rmse=test_rmse, r2=test_r2),
        )
    )
end
