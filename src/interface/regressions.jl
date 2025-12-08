function fit_growth_models(
    groupings::OrderedDict{String,DataFrame}; degree::Int=2
)::PolyGrowthModel
    models = PolyGrowthFunction[]

    g = length(values(groupings))
    m = length(ALL_METRICS)

    train_ = NamedTuple((
        Symbol(k) => Vector{Float32}(undef, g) for k in CoralFlow.ALL_METRICS
    ))
    test_ = deepcopy(train_)

    for (i, (_, df)) in enumerate(groupings)
        isempty(df) && continue

        sub_df = df[df[!, TRAIN_CLASS] .> 0, :]

        x_idx = sortperm(sub_df.diam)
        xi = sub_df.diam[x_idx]
        yi = sub_df.diamnext[x_idx]

        # Fit model
        m = curve_fit(Polynomial, log.(xi), log.(yi), degree)
        # m = curve_fit(Polynomial, xi, yi, degree)

        # Create growth function
        model = PolyGrowthFunction(xi, yi, m)
        prediction = model.(xi)
        push!(models, model)

        # Collate training metrics
        for m in ALL_METRICS
            getfield(train_, Symbol(m))[i] = m(prediction, yi)
        end

        # As above, but for test data
        sub_df = df[df[!, TEST_CLASS] .> 0, :]
        x_idx = sortperm(sub_df.diam)
        xi = sub_df.diam[x_idx]
        yi = sub_df.diamnext[x_idx]

        prediction = model.(xi)

        # Collate training metrics
        for m in ALL_METRICS
            getfield(test_, Symbol(m))[i] = m(prediction, yi)
        end
    end

    return PolyGrowthModel(
        collect(keys(groupings)),
        models,
        (;
            train=train_,
            test=test_
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
function fit_survival_models(
    groupings::OrderedDict{String,DataFrame}; degree::Int64=2
)::PolySurvivalModel
    models = PolySurvivalFunction[]
    g = length(values(groupings))
    m = length(ALL_METRICS)

    train_ = NamedTuple((
        Symbol(k) => Vector{Float32}(undef, g) for k in CoralFlow.ALL_METRICS
    ))
    test_ = deepcopy(train_)

    for (i, (_, df)) in enumerate(groupings)
        isempty(df) && continue

        sub_df = df[df[!, TRAIN_CLASS] .> 0, :]

        # Fit model based on observed size at time of mortality (`diam_mort`)
        x_idx = sortperm(sub_df.diam_mort)
        xi = sub_df.diam_mort[x_idx]
        yi = sub_df[x_idx, TRAIN_CLASS_MEAN_ID]

        # Fit model
        m = curve_fit(Polynomial, log.(xi), yi, degree)

        # Create survival function
        model = PolySurvivalFunction(xi, yi, m)

        # Calculate metrics
        prediction = model.(xi)
        push!(models, model)

        # Collate training metrics
        for m in ALL_METRICS
            getfield(train_, Symbol(m))[i] = m(prediction, yi)
        end

        # Repeat above for test data
        sub_df = df[df.class_test .> 0, :]
        x_idx = sortperm(sub_df.diam_mort)
        xi = sub_df.diam_mort[x_idx]
        yi = sub_df.class_test_mean[x_idx]
        prediction = model.(xi)

        # Collate training metrics
        for m in ALL_METRICS
            getfield(test_, Symbol(m))[i] = m(prediction, yi)
        end
    end

    return PolySurvivalModel(
        collect(keys(groupings)),
        models,
        (;
            train=train_,
            test=test_
        )
    )
end
