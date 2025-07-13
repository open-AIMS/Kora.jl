RMSE(y_hat, y) = sqrt(mean((y_hat .- y) .^ 2))
R2(y_hat, y) = 1 - (sum((y .- y_hat) .^ 2) / sum((y .- mean(y)) .^ 2))

pearson(y_hat, y) = StatsBase.cor(y_hat, y)
spearman(y_hat, y) = StatsBase.corspearman(y_hat, y)
kendall(y_hat, y) = StatsBase.corkendall(y_hat, y)

const ALL_METRICS = [RMSE, R2, pearson, spearman, kendall]
